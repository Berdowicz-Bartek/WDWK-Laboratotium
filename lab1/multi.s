EXIT_NR 	= 1
READ_NR 	= 3
WRITE_NR 	= 4

STDOUT			= 1
EXIT_CODE_SUCCESS	= 0

.data
	n:	.int 10
.text
	
.global _start

_start:

	movl n, %eax
	imull %eax, %eax
	movl %eax, %ebx
	imull n, %eax
	movl %eax, %ecx		#ecx = n^3
	
	movl %ebx, %eax
	imull $3, %eax	#eax = 3n^2
	
	addl %eax, %ecx
	
	movl n, %eax
	imull $2, %eax 		#eax = 2n
	
	subl %eax, %ecx		#ecx = n^3 + 3n^2 - 2n
	
	

	mov $EXIT_NR, %eax
	mov $EXIT_CODE_SUCCESS, %ebx
	int $0x80
