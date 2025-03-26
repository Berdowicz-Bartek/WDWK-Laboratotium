.section .note.GNU-stack,"",%progbits

.text
    .global foo
    .extern printf

foo:
    rdtsc
    mov %rdx, %r8
    mov %rax, %r9

    rdtsc
    sub %r8, %rdx
    sbb %r9, %rax
    
    mov $fmt, %rdi
    mov %r9, %rsi
    xor %rax, %rax
    call printf

    ret

.section .rodata
fmt:
    .asciz "Czas wykonania: %llu cykli\n"
