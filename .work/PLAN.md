# Shortest-Path Competition 優化計畫（2026-09-19）

給執行模型的說明：這份文件是完整的工作指示。每個工作包（WP）都可以獨立交給一個代理，
共用的規則、資料、指令在 `.work/BRIEF.md`，開工前先讀它，再讀整份 `solver.cpp`。
所有量測都在 WSL 進行（`wsl -e bash -lc '...'`），資料在 `/root/spc/instances/`。

---

## 0. 目標與底線

- 目標：把 PYPARTY 的**真實** overall geomean 推到最高；排行榜第一名目前是 442x。
- 現況（安靜機器、grade.py 真實量測，檔案 `/root/spc/result.json`）：

| 實例 | T_base (s) | 現在 T_solver (s) | 加速 | 分段 (s) |
|---|---:|---:|---:|---|
| road2d    | 1013 | 15.5 | 65   | contract 9.6 / hub labels 6.3 / queries 0.15 |
| lattice3d | 1526 | 30.2 | 51   | contract ~25–35 / CH queries 8 |
| local2d   | 596  | 10.0 | 60   | contract 4.9 / CH queries 9.8（33 µs/筆） |
| scalefree | 1378 | 15.3 | 90   | bidir Dijkstra queries 全部 |
| hugeq     | 1764 | 1.3  | 1359 | contract 1.28 / labels 0.58 / queries 0.42 |
| wide64    | 882  | 41.5 | 21   | contract ~35 / CH queries ~10 |
| **overall** | | | **89.5** | |

- 要贏 442x，六個實例的乘積要達到 442^6 ≈ 7.5e15。若某實例只能到 30x，其餘五個平均要 790x。
  這在單執行緒下極難；計畫的定位是「合法手段能拿多少拿多少」，每翻一倍就等值於課程曲線上的一步。
- 硬性規則（README 規則 1–8）：單一 `solver.cpp`、C++17、只准標準函式庫標頭、單執行緒（grade.py
  以 cpu > 1.4×wall 判定 THREADS）、<4 GB、<3600 s、輸出與 foundation 完全一致、不得內含預先算好的答案。
  可用：`#pragma GCC optimize/target`、GCC vector extension、inline asm（不需要額外標頭）。
  不可用：GPU、多執行緒、`<immintrin.h>`、`<sys/mman.h>`、任何依賴發布檔案「確切位元組」的東西
  （最終評分用同參數、不同 seed 的新實例）。
- **絕對禁止**：修改 `result.json` 的任何數字、修改 `grade.py`、在跑 baseline 時讓機器變忙。
  這些是資料造假；教授手上有 solver.cpp，重跑一次就會露餡。只交 grade.py 自己寫出的 result.json。

---

## 1. 為什麼是這幾個方向（瓶頸分析）

1. **收縮（CH contraction）是最大成本**：road2d 9.6 s、lattice3d ~30 s、wide64 ~35 s、local2d 4.9 s。
   現在用「懶惰優先佇列 + simulate() 估算 + 有界見證搜尋」，在網格（lattice）上特別差，因為
   網格沒有天然的階層，懶惰排序會產生大量捷徑與見證搜尋（wide64 見證搜尋 settle 2 億次）。
2. **網格類查詢太慢**：local2d 每筆 33 µs、lattice3d 每筆 270 µs。原因同上：排序差 → 向上搜尋空間大。
3. **scalefree 完全靠雙向 Dijkstra**：每筆 300 µs。power-law 圖上幾乎所有最短路徑都經過少數 hub，
   PLL（pruned landmark labeling）在這類圖上標籤只有幾十筆，查詢 <1 µs。
4. **I/O 不是瓶頸**：`oracle_solver.cpp`（只讀檔、抄答案）在 dev tier 2–7 ms；大檔案 (100 MB) 解析約 0.1 s。

---

## 2. 工作包（可四路平行；每包一個代理、一個檔案、一個 WSL 目錄）

### WP-lattice：巢狀分割順序的 CH（檔案 `.work/lat.cpp`，TAG=lat）
對象：lattice3d、local2d、wide64（無向網格）。目標：lattice3d contract <4 s、total <8 s；
local2d contract <0.8 s；wide64 contract <8 s。

步驟：
1. **偵測網格**：解析後檢查是否為 d 維環面（torus）網格且編號為 row-major：找 side 使
   V = side^d（d=2 或 3），每個頂點 i 的鄰居恰為 i±1（同列、環繞）、i±side、i±side²。
   偵測失敗就走原本的懶惰排序（road2d/hugeq 必須不受影響）。
2. **產生 nested-dissection 順序**：對（子）網格遞迴二分，永遠切最長軸；環面的軸在最頂層要切
   兩片平行超平面才能真正分開，之後就是普通的盒子。分隔面頂點拿該區域最高的 rank，兩邊子區域
   遞迴拿較低的 rank；葉子（≤32–64 頂點）順序任意或用簡單 edge-difference。輸出 `rank_of[v]`。
3. **依固定順序收縮**：沒有優先佇列、沒有 simulate()、沒有懶惰更新；每個頂點依 rank 順序呼叫
   find_shortcuts（有界見證搜尋，調 CON_SETTLE，可加 hop 上限），從動態圖移除、加捷徑，與現有迴圈相同。
   評估拿掉 `prune_redundant_arcs`。
4. 用現有 `ch_query` 量查詢階段；回報 upward arc 數與每筆查詢時間。
5. 若有餘裕：縮小見證搜尋成本（stamp 重用、小 heap）、葉子用更便宜的順序。

驗收：`devcheck.sh lat.cpp` 六個 OK；`bench.sh` 三個大實例 OK 且分段時間達標或明確回報差距。

### WP-local：短程查詢快速路徑（檔案 `.work/loc.cpp`，TAG=loc）
對象：local2d（90% 查詢的 Dijkstra rank 在 2^4..2^12）、wide64（rank 2^8..2^15）。
目標：兩者查詢階段皆 <1.5 s。只改查詢相關的新函式與 run_queries/main 的分派，不碰收縮程式碼
（會被合併進 WP-lattice 的檔案）。

步驟：
1. 另建一份原圖的 CSR（鄰接依權重排序），與 CH 並存（記憶體充裕）。
2. **有 settle 預算的雙向 Dijkstra**：預算 N（約 2000–8000，可調）內兩邊相遇且 min-key 和 ≥ mu 就回傳，
   否則退回 `ch_query`。用 timestamp 代替清陣列、重用 heap、查詢內零配置。
   local2d 權重 ≤1000 → 用 Dial 桶佇列；wide64 64 位元距離 → radix heap。
3. 依量測調預算；回報快速路徑命中率與查詢階段時間。
4. 優化剩下的 CH 查詢：stall 檢查與鬆弛合併成一次掃描、用 mu 剪枝、更便宜的 heap。
5. 選做：所有查詢事先已知，試按來源編號（即網格 row-major）排序處理以改善快取，輸出仍按原順序。

驗收：`devcheck.sh loc.cpp` 六個 OK；local2d、wide64 大實例 OK；回報新增了哪些函式/行號供合併。

### WP-road：道路圖收縮加速與 HL/CH 抉擇（檔案 `.work/road.cpp`，TAG=road）
對象：road2d（20 萬筆 rank 2^6..2^15 查詢）、hugeq（60 萬筆均勻查詢、小圖）。
目標：road2d total <3 s、hugeq total <0.6 s。收縮程式碼也被其他家族共用，收縮變快全體受益。

步驟：
1. 剖析收縮：`prune_redundant_arcs`、初始優先權（每個頂點跑 simulate）、懶惰重估、見證搜尋
   （g_searches/g_settles 有印出）、動態圖 `vector<vector<DArc>>` 的線性掃描各占多少。
   試 `PRUNE_SETTLE=0`、更低的 CON_SETTLE/SIM_SETTLE、hop 上限；見證搜尋只需證明路徑 ≤ 捷徑權重，
   所有目標解決或 key 超過最大待決門檻就停。
2. RoutingKit 式排序：優先權 = edge difference + 已收縮鄰居數（level）+ hop quotient，只重估被收縮
   頂點的鄰居，simulate 更便宜。確認 upward arc 數與查詢時間沒變差（道路捷徑 ≈1.2×E 正常）。
3. 以量測決定每個實例用 hub labels 還是 CH 查詢：road2d 查詢是短程，20 萬筆雙向 CH 查詢可能 2–3 s，
   而 labels 建置 6.3 s；hugeq 60 萬筆查詢在小圖上 labels 明顯划算。抉擇寫成 (V,E,Q) 的啟發式，
   新 seed 的實例也要成立。可先試更快的 label 建置（`HL_PRUNE` 剪枝成本）。
4. 微優化：CSR 式動態圖、避免重複 INF 檢查、距離放得下就用 32 位元、解析 37 MB 檔案時整段跳過座標區塊。

驗收：`devcheck.sh road.cpp` 六個 OK；road2d、hugeq 大實例 OK；回報 HL-vs-CH 規則與各旋鈕最佳值。

### WP-scalefree：PLL 與 hub 目標批次化（檔案 `.work/sf.cpp`，TAG=sf）
對象：scalefree（有向 power-law，3% 為 sink 所以約 3% 查詢答 -1；70% 均勻、30% 目標是 hub）。
目標：total <4 s（=345x）。

步驟：
1. 先量現有未啟用的 PLL 路徑：`USE_PLL=1 USE_BIDIR=0`（build_pll/pll_query），先跑 scalefree_dev。
   記錄建置時間、標籤大小、查詢時間。建置太慢就優化：root 依 (in+1)*(out+1) 遞減、CSR 鄰接、
   當前距離已被既有標籤回答時停止擴展、32 位元距離、timestamp、有向圖對 sink/source 的處理。
2. 獨立改進雙向 Dijkstra 作為後備：把 30% hub 目標查詢**按目標分組**，每個目標只跑一次反向 Dijkstra，
   該目標的所有來源都 settle 就停；均勻查詢維持雙向。檢查平衡規則與依權重排序的提前中斷是否浪費工作。
3. 以量測選方法，規則為「有向且無座標 → 此路徑」。
4. 微優化：14.6 MB 檔案解析、記憶體配置、radix/4-ary heap、避免重複清陣列。

驗收：`devcheck.sh sf.cpp` 六個 OK；scalefree 大實例 OK（含 -1）；回報標籤統計與選擇規則。

---

## 3. 整合與交付（由整合者執行，順序固定）

1. 合併順序：以 `lat.cpp`（收縮改動最大）為底 → 併入 `road.cpp` 的排序/HL 改動 → 併入 `loc.cpp`
   的查詢函式與分派 → 併入 `sf.cpp` 的 PLL/批次路徑。每併一步跑一次 `devcheck.sh`。
2. 在合併後檔案頂端加：`#pragma GCC optimize("O3")` 與 `#pragma GCC target("avx2,bmi,bmi2,popcnt")`
   （教授端編譯旗標未知，寫進原始碼才能在最終評分也生效）。用 devcheck 確認沒有壞掉。
3. 合併檔複製為根目錄 `solver.cpp`，關掉其他程式，執行 `bash run.sh large`（約 30–60 分鐘），
   讀回真實的 `result.json`。任何實例 WRONG/NONDETERMINISTIC 都要先修，不得帶病上傳。
4. 用**同一次跑出的** `solver.cpp` + `result.json` 上傳到記分板（sha256 要對得上）。
5. 更新 `.work/BRIEF.md` 的現況表與記憶檔，記錄新的分段時間，方便下一輪。

---

## 4. 若時間還有剩：第二輪候選（依預期收益排序）

1. lattice3d 查詢：ND 順序下用 CH 的「many-to-many 桶」或查詢剪枝，目標每筆 <100 µs。
2. wide64：偵測到網格後可推出座標，用 `w_min × 網格距離` 當 A* 下界（權重 1e8..1e9，下界約真值 1/3），
   對 rank 查詢有幫助。
3. road2d：座標存在時用 ALT（landmark）或 CH+A*；20 萬筆短程查詢可能比 HL 更划算。
4. 解析：一次 fread 後以 8 位元組一組的 SWAR 解析數字；輸出用查表法轉字串。
5. hugeq：60 萬筆查詢按來源分組，label 一次載入快取。

---

## 5. 風險

- 網格偵測或 ND 切法對「同參數、不同 seed」的新實例必須同樣成立；side 由 V 推算，不要寫死。
- 任何方法都要保留退路（偵測失敗就走舊路徑），避免最終評分出現 WRONG → 0.1x。
- 四個代理同時跑 benchmark，時間噪音 ±15%；最終計分一定要在安靜的機器上獨跑。
- 記憶體：wide64 現在 ~660 MB；加 CSR/標籤後仍需 <4 GB（grade.py 用 RLIMIT_AS 強制）。

---

## 6. 測資特化：對這六個家族特別快、對任何輸入仍然正確

原則：**在載入時偵測結構 → 通過驗證才走特化路徑 → 否則走通用路徑**。所有特化路徑都是精確演算法
（輸出必須與 foundation 逐位元相同），偵測失敗只會變慢不會變錯。README 明說「依賴圖結構的手法會
延續到新 seed 的實例」，所以這些都合法；唯一禁止的是依賴發布檔案的確切內容。

### 6.1 結構偵測（載入後、任何演算法前，總成本 <50 ms）

| 偵測項目 | 判定方式 | 命中的實例 | 解鎖的特化 |
|---|---|---|---|
| d 維環面網格、row-major 編號 | 由 V 推 side（V = side^2 或 side^3），逐邊驗證每條邊都是 i↔i±1（同列環繞）/ i±side / i±side² 且每個頂點恰好 2d 個鄰居 | local2d、wide64、lattice3d | 隱式網格、ND 順序、座標下界 |
| 網格狀道路圖（有座標） | 座標區塊存在；不需驗證，只把座標當啟發式用 | road2d、hugeq | ND 順序（由座標二分）、座標下界 |
| 小權重 | 解析時記下 w_max；w_max ≤ 約 4096 | local2d | Dial 桶佇列 |
| 32 位元距離 | 由 (V, w_max) 或 meta 推不出來時，維持 64 位元；只有 w_max×直徑上界 < 2^32 時才用 32 位元 | road2d、hugeq、local2d、scalefree、lattice3d | 距離陣列減半、快取更好 |
| 有向、無座標、度數 power-law | 已有的 USE_BIDIR 條件，再加「最大入度 > 100×平均度」 | scalefree | PLL / hub 特化 |
| 查詢集統計 | 讀完查詢檔就有：Q、相異來源數、相異目標數、s==t 數、重複 (s,t) 數 | 全部 | 6.3 節的查詢集特化 |

### 6.2 各家族的特化路徑（含退路與預期收益）

**A. 隱式網格（local2d、wide64、lattice3d）**
- 不建 CSR：鄰居用算術得到（i±1 含環繞、i±side、i±side²），只存每個頂點往 +x/+y/+z 的權重
  （d 個 `uint32/uint64` 陣列）。省掉索引載入與一半記憶體，Dijkstra 每次 settle 的成本降到約 20–40 ns。
- ND 順序直接由 (x,y,z) 算，不需要任何圖搜尋（WP-lattice 第 2 步的輸入）。
- 退路：偵測失敗 → 原本的 CSR 與懶惰排序。

**B. wide64：不建 CH，直接答（最大單項收益，41.5 s → 估 4–8 s）**
- 查詢全是 rank 2^8..2^15 的短程查詢，而收縮要 35 s。做法：對每筆查詢跑「有 settle 預算」的雙向
  Dijkstra（隱式網格 + radix heap，64 位元 key）。**CH 改為惰性建置**：只有某筆查詢超過預算
  （例如 2^17 settles）才建 CH 並改用它答剩下的查詢。在這個家族上 CH 永遠不會被建，但任何
  其他輸入仍然正確。
- 加速下界：偵測到網格後有座標，用 `w_min × 環面曼哈頓距離` 當 A* 下界（w_min 由解析得到，
  loguniform 1e8..1e9 的 w_min 約是平均權重的 1/4，能省 20–30% settles）。雙向 A* 用平均勢能
  （(h_f − h_b)/2）保持一致性。
- 再進一步：ALT（8 個 landmark，隱式網格上每個 Dijkstra 約 0.15 s，共 1.2 s 前處理），下界緊很多，
  短程查詢的 settles 可再降 3–5 倍。估總時間 4 s（約 220x）；純雙向 Dijkstra 估 7–9 s（約 100–125x）。
- 退路：預算超過 → 惰性建 CH；網格偵測失敗 → 通用路徑。

**C. local2d：Dial 桶 + 惰性 CH**
- 90% 查詢 rank ≤ 4096、權重 ≤ 1000：隱式網格 + Dial 桶佇列（環狀 1024 桶，push/pop O(1)）的雙向
  搜尋，估每筆 5–10 µs → 27 萬筆約 1.5–2.5 s。
- 10% 均勻查詢（3 萬筆）在 30 萬頂點網格上雙向 Dijkstra 要約 3 ms/筆，太慢，所以 ND-CH 仍要建
  （隱式網格上估 0.5–0.8 s），只用來答超過預算的查詢。
- 預期：解析 0.05 + ND-CH 0.7 + 本地 2 + 全域 0.6 ≈ 3.3 s（約 180x）。

**D. lattice3d：ND-CH 為主，ALT 為備援**
- 30 萬頂點 3D 環面、3 萬筆全域查詢。主線是 WP-lattice 的 ND-CH；3D 的頂層分隔面有 67² ≈ 4500 個頂點，
  其間捷徑可能到數百萬條，要量。若收縮或查詢爆掉，備援是隱式網格 + 16 個 landmark 的雙向 ALT
  （前處理約 0.5 s，估每筆 0.5–1 ms → 15–30 s，仍比現在的 30–45 s 好）。

**E. road2d / hugeq：座標 ND 順序 + 查詢集分組**
- 生成器把道路圖建在 548²/224² 的網格上（arterial stride 8、highway stride 64、8% 單向），編號同樣
  row-major。可直接用座標做 ND 二分當收縮順序，跳過整個優先佇列與 simulate()（收縮估 9.6 s → 1–2 s）；
  與重要度排序相比查詢空間可能稍大，用 WP-road 第 3 步的量測決定 HL 或 CH。
- hugeq：Q=60 萬、V=5 萬，每個來源平均出現 12 次。把查詢依來源分組，同一來源的 forward label
  只載入一次（HL）或 forward 向上搜尋只跑一次（CH one-to-many）。
- 退路：無座標 → 現有懶惰排序。

**F. scalefree：PLL 為主，hub 上界為備援**
- 主線：WP-scalefree 的 PLL。標籤依 hub rank 排序，查詢是兩個短列表的合併（可向量化）。
- 備援（PLL 建置太慢時）：取入度最高的 K=16 個 hub，各跑一次正向與反向 Dijkstra（約 1 s），
  每筆查詢先用 `min_h d(s,h)+d(h,t)` 得上界 U，再跑以 U 剪枝的雙向 Dijkstra；因為真實最短路徑幾乎
  都經過 hub，U 通常就是答案，搜尋很快收斂。任何圖上都正確。
- 3% sink：來源是 sink 且 s≠t → 直接 -1（已存在），目標入度 0 → -1。

### 6.3 查詢集特化（六個實例通用，零風險）
- 讀完整個查詢檔後再開始計算（現在已是如此），先做：s==t → 0；完全相同的 (s,t) 只算一次；
  依來源排序處理以提升快取命中，輸出時還原原順序。
- 「先用小預算的便宜方法，失敗再用重方法」的兩段式查詢（B、C 的核心），讓前處理可以惰性化：
  **只有真的需要時才付前處理的錢**。
- 解析特化：檔案只含非負整數與空白，用 8 位元組一組的 SWAR 找分隔符與轉整數；道路圖的座標區塊
  只在需要 ND/座標下界時解析，否則整段跳過。

### 6.4 納入工作包的方式
- WP-lattice 增加：6.1 的網格偵測與 6.2-A 的隱式網格（其他兩包都靠它）。
- WP-local 增加：6.2-B 的惰性 CH 與 A*/ALT 下界、6.2-C 的 Dial 桶；量測「純雙向 / A* / ALT」三種。
- WP-road 增加：6.2-E 的座標 ND 順序、hugeq 的來源分組。
- WP-scalefree 增加：6.2-F 的 hub 上界備援。
- 整合階段增加：6.3 的查詢去重、排序、SWAR 解析。

### 6.5 加入特化後的預期

| 實例 | 第 2 節樂觀 | 加入第 6 節 | 主要來源 |
|---|---:|---:|---|
| road2d | 338x (3 s) | 400x (2.5 s) | 座標 ND 收縮 |
| lattice3d | 191x (8 s) | 191x (8 s) | 不變（ND-CH 已是主線） |
| local2d | 248x (2.4 s) | 180–250x (2.4–3.3 s) | Dial + 惰性 CH，看本地查詢實測 |
| scalefree | 345x (4 s) | 345x (4 s) | 不變 |
| hugeq | 2940x (0.6 s) | 3500x (0.5 s) | 來源分組 |
| wide64 | 88x (10 s) | 150–220x (4–6 s) | 不建 CH + ALT 下界 |
| **overall** | **335x** | **約 380–430x** | |

仍略低於 442x；能不能過取決於 wide64 的 ALT 與 local2d 的本地查詢實測，這兩項的估計誤差最大。
