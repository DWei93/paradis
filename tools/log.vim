" Vim syntax file
" Language: *.log syntax
" Maintainer: Dean Wei
" Latest Revision: 6 February 2017

if exists("b:current_syntax")
   finish
endif

syn keyword item Last git log Modification Running Info Job Message
syn keyword highlight1 '@@@' '-->'
syn keyword highlight2 Failed
syn keyword highlight3 SUCCESS

"syn match number '\d\+'
"syn match number '\d\+'
"syn match number '[-+]\d\+'

let b:current_syntax = "cel"

hi def link number        Constant
hi def link highlight1    Comment
hi def link highlight2    Type 
hi def link highlight3    PreProc
hi def link item          Statement

