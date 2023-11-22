void main()
{
    char *var = "foo:bar";
    char *subst = "quux";
    char *new;
    str_parse(var, "s/^(.+?):(.+)$/$1-%s-$2/", &new, subst);
    /* now we have: var = "foo:bar", new = "foo:quux:bar" */
    printf("'foo:quux:bar' ?= '%s'\n",new);
    free(new);
    return;
}
