proc iml;
start MyFunc(x);
   return( x );
finish;
/** integrate cos(x) on [a,b] **/
a = {0, 1, 2, 3};    
b = {3, 0, 5, 4};  
call quad(R, "MyFunc", a||b);
/** compare with exact answer */
exactInt = sin(b) - sin(a);
diff = exactInt - R;
print R exactInt diff;
create test var {"R" "a" "b"};
append;
close test;
