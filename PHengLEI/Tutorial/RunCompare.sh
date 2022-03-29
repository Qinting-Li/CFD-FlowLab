#cd PHengLEIv3d0
#cd Tutorial

ProjectName=$3
ExeName=$4
ExeOutputFile=$5
StandardOriginFile=$6

cd $3
while [ -n "$1" ]
do
	case $1 in
	-a)
		pwd
		echo "input file1 parameter is : $ProjectName"
		echo "input file2 parameter is : $ExeName"
		echo "input file2 parameter is : $ExeOutputFile"
		echo "input file3 parameter is : $StandardOriginFile"
		pwd
		if [ -e "in" ];then #-e 表示判断文件夹是否存在，不存在运行else后面的部分。
			if [ -e "out" ];then
				#echo "folder in and out exist!"
				echo "$ExeName"
				$ExeName 
				../CompareFile.sh -a -- $ExeOutputFile $StandardOriginFile 
			else
				echo "file out not exist!"
				cd ../ 
				
				pad=$(printf '%0.1s' "-"{1..49})
				padlength=100
				string1="folder out not exist!"
				string2='Skip'
				printf "\n"
				printf '%s' "$string1"
				printf '%*.*s' 0 $((padlength - ${#string1} - ${string2} )) "$pad"
				printf "\033[37;33m $string2 \033[0m"
				printf "\n"
		
			fi
		else
			#echo "folder in not exist!"
			cd ../ 
		
			pad=$(printf '%0.1s' "-"{1..49})
			padlength=100
			string1="folder in not exist!"
			string2='Skip'
			printf "\n"
			printf '%s' "$string1"
			printf '%*.*s' 0 $((padlength - ${#string1} - ${string2} )) "$pad"
			printf "\033[37;33m $string2 \033[0m"
			printf "\n"
		fi
	;;
	
	-s)
		if [ -e "in" ];then #-e 表示判断文件夹是否存在，不存在运行else后面的部分。
			if [ -e "out" ];then
				$ExeName > /dev/null
				../CompareFile.sh -s -- $ExeOutputFile $StandardOriginFile 
			else
				echo "file out not exist!"
				cd ../ 
				
				pad=$(printf '%0.1s' "-"{1..49})
				padlength=100
				string1="folder out not exist!"
				string2='Skip'
				printf "\n"
				printf '%s' "$string1"
				printf '%*.*s' 0 $((padlength - ${#string1} - ${string2} )) "$pad"
				printf "\033[37;33m $string2 \033[0m"
				printf "\n"
		
			fi
		else
			#echo "folder in not exist!"
			cd ../ 
		
			pad=$(printf '%0.1s' "-"{1..49})
			padlength=100
			string1="folder in not exist!"
			string2='Skip'
			printf "\n"
			printf '%s' "$string1"
			printf '%*.*s' 0 $((padlength - ${#string1} - ${string2} )) "$pad"
			printf "\033[37;33m $string2 \033[0m"
			printf "\n"
		fi
	;;
	
	--) shift
		break;;
		
	*) echo "&&&&&&&&&&&&&&&&";;
	esac
	shift
done
#pwd

