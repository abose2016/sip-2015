#include <stdio.h>
int addTwoNumbers(int, int);
main()
{
 int x = 0;
 int y = 0;
Chapter 7 • Pointers 159
 printf("\nEnter first number: ");
 scanf("%d", &x);
 printf("\nEnter second number: ");
 scanf("%d", &y);
 printf("\nResult is %d\n", addTwoNumbers(x, y));
} //end main
int addTwoNumbers(int x, int y)
{
 return x + y;
} //end addTwoNumbers
