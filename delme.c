#include "tfm.h"

int main(void)
{
    fp_int a;
    char buf[4096];

    fp_init(&a);
    fp_read_radix( &a,
        "///////////93zgY8MZ2DCJ6Oek0t1pHAG9E28fdp7G22xwcEnER8b5A27cED0JT"
        "xvKPiyqwGnimAmfjybyKDq/XDMrjKS95v8MrTc9UViRqJ4BffZVjQml/NBRq1hVj"
        "xZXh+rg9dwMkdoGHV4iVvaaePb7iv5izmW1ykA5ZlmMOsaWs75NJccaMFwZz9CzV"
        "WsLT8zoZhPOSOlDM88LIkvxLAGTmbfPjPmmrJagyc0JnT6m8oXWXV3AGNaOkDiux"
        "uvvtB1WEXWER9uEYx0UYZxN5NV1lJ5B9tYlBzfLO5nWvbKbywfLgvHNI9XYO+WKG"
        "5NAEMeggn2sjCnSD151wCwXL8QlV7BfaxFk515ZRxmgAwd5NNGOCVREN3uMcuUJ7"
        "g/MkZDi9CzSUZ9JWIYLXdSxZqYOQqkvhyI/w1jcA26JOTW9pFiXgP58VAnWNUo0C"
        "k+4NLtfXNMnt2OZ0kjb6uWZYJw1qvQinGzjR/E3z48vBWj4WgJhIol//////////",
        64 );

    if( fp_isprime( &a ) ) printf("It's prime.\n");
    else printf( "Not prime.\n");

    return 0;
}
