
void
aprint_
( char * sz , int * len )
{
    int l;
    for (l=0; l<*len; l++)
    {
        printf(" %d",*(sz+l));
    }
    printf("\n");
    return;
}

