function PrintfLine_same {
		pad=$(printf '%0.1s' "-"{1..49})
		padlength=100
		string2='true'
		string10=$file1
		string11=' is the same thing as '
		string12=$file2
		string1=${string10}${string11}${string12}
		printf "\n"
		printf '%s' "$string1"
		printf '%*.*s' 0 $((padlength - ${#string1} - ${string2} )) "$pad"
		printf "\033[32;32m true \033[0m"
		printf "\n"
}

function PrintfLine_diff {
		pad=$(printf '%0.1s' "-"{1..49})
        padlength=100
        string2='false'
        string10=$file1
        string11=' is not the same thing as '
        string12=$file2
        string1=${string10}${string11}${string12}

		printf "\n"
        printf '%s' "$string1"
        printf '%*.*s' 0 $((padlength - ${#string1} - ${string2} )) "$pad"
		printf "\033[31;31m false \033[0m"
		printf "\n"
}
function PrintfLine_NoExist {
	pad=$(printf '%0.1s' "-"{1..49})
	padlength=100
	string2='Skip'
	string1="$file1 or $file2 does not exitst,please check filename."
	printf "\n"
	printf '%s' "$string1"
	printf '%*.*s' 0 $((padlength - ${#string1} - ${string2} )) "$pad"
	printf "\033[32;36m $string2 \033[0m"
	printf "\n"
}

 
while [ -n "$1" ]
do
	case $1 in
	-a)
		echo "Run command: -$0 file1 file2"
		
		file1=$3
		file2=$4
		
		echo "input file1 name is : $file1"
		echo "input file2 name is : $file2"
		
		if [ -f $file1 ] && [ -f $file2 ]
		then 
			diff -w $file1 $file2 > /dev/null
			if [ $? != 0 ]
			then
				PrintfLine_diff
			else
				PrintfLine_same
			fi
		else
			# $file1 or $file2 does not exitst.
			PrintfLine_NoExist
		fi
	;;
	
	-s)
		file1=$3
		file2=$4
		if [ -f $file1 ] && [ -f $file2 ]
		then 
			diff -w $file1 $file2 > /dev/null
			if [ $? != 0 ]
			then
				PrintfLine_diff
			else
				PrintfLine_same
			fi
		else
			# $file1 or $file2 does not exitst.
			PrintfLine_NoExist
		fi
	;;
	
	--) shift
		break
	;;
	
	*) echo "check input option"
	;;
	
	esac
	shift
done
