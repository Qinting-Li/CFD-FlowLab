import platform
import os

if platform.system().lower() == 'windows':
	print("win")
	cmd = "Run_UtilTest_windows.sh -a"
	val = os.system(cmd)
	print(val)
elif platform.system().lower() == 'linux':
	print("lin")
	cmd = "Run_UtilTest_linux.sh"
	val = os.system(cmd)
	print(val)
