1.
	ExecuteShellScript.py 负责区分Linux和windows环境，并使用指令执行shell脚本。
	Run_UtilTest_linux.sh 在linux下负责 找到并按照指令执行CompareFile.sh脚本
	Run_UtilTest_windows.sh 在windows下负责 找到并按照指令执行CompareFile.sh脚本
	CompareFile.sh 负责比对 单元测试输出结果 和 标准结果 是否相同，相同就输出绿色的true，不同就输出红色的false。同时打印出是什么文件在进行比对。

2.
	在ExecuteShellScript.py的cmd参数中，可以使用两种指令。以windows下为例，
	第一种是全输出，Run_UtilTest_windows.sh -a。
	第二种是简要输出，Run_UtilTest_windows.sh -s。
	对应linux一样。

