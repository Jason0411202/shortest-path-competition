s = open('sp.cpp', 'rb').read().decode()
s = s.replace('%llu touched\n", b,', '%llu touched\\n", b,')
open('sp.cpp', 'wb').write(s.encode())
