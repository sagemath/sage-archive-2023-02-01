
#
# SAGE support utilities to read into the GAP session.
#
\$SAGE := rec();

\$SAGE.OldPager := Pager;


\$SAGE.NewPager :=
          function( data )
    local   str,  lines,  line, fn, start;
    str := OutputTextFile(\$SAGE.tempfile,false);
    start := 1;
    if IsRecord(data) then
        lines := data.lines;
        if IsBound(data.start) then
            start := data.start;
        fi;
    else
        lines := data;
    fi;
    if IsString(lines) then
        lines := SplitString(lines,"\n");
    fi;
    for line in lines do
        WriteLine(str, line);
    od;
    CloseStream(str);
    Print("Page from ",start,"\n");
end;

\$SAGE.StartInteract := function()
    MakeReadWriteGlobal("Pager");
    Pager := \$SAGE.OldPager;
    HELP_VIEWER_INFO.screen.show := \$SAGE.OldPager;
    MakeReadOnlyGlobal("Pager");
end;


\$SAGE.StopInteract := function()
    MakeReadWriteGlobal("Pager");
    Pager := \$SAGE.NewPager;
    HELP_VIEWER_INFO.screen.show := \$SAGE.NewPager;
    MakeReadOnlyGlobal("Pager");
end;


\$SAGE.StopInteract();

#\$SAGE.ErrorHandler := function(m,a,m2,mode)
#    PrintTo("*errout*", m);
#    if a <> fail then
#        PrintTo("*errout*",a);
#    fi;
#    SetErrorHandler(\$SAGE.ErrorHandler);
#    return true;
#end;

#SetErrorHandler(\$SAGE.ErrorHandler);

SetAllInfoLevels(0);

\$SAGE.OperationsAdmittingFirstArgument := function(obj)
    local   hits,  myflags,  i,  flagss,  flags;
    hits := [];
    myflags := FlagsType(TypeObj(obj));
    for i in [1,3..Length(OPERATIONS)-1] do
        flagss := OPERATIONS[i+1];
        for flags in flagss do
            if Length(flags) >= 1 and IS_SUBSET_FLAGS(myflags, flags[1]) then
                Add(hits, OPERATIONS[i]);
                break;
            fi;
        od;
    od;
    return hits;
end;


\$SAGE.CleanOperationName := function(name)
    local   lt,  ls;
    lt := Length("Tester(");
    if Length(name) > lt and name{[1..lt]} = "Tester(" then
        return Concatenation("Has",name{[lt+1..Length(name)-1]});
    fi;
    ls := Length("Setter(");
    if Length(name) > ls and name{[1..ls]} = "Setter(" then
        return Concatenation("Set",name{[lt+1..Length(name)-1]});
    fi;
    return name;
end;

\$SAGE.HasAtLeastOneMethodAsFirstArgument := function(op,obj)
    local   t,  f,  n,  meths,  i;
    t := TypeObj(obj);
    f := FlagsType(t);
    for n in [1..6] do
        meths := METHODS_OPERATION(op,n);
        for i in [1,n+5..LENGTH(meths)-n-3] do
            if IS_SUBSET_FLAGS(f,meths[i+1]) then
                return true;
            fi;
        od;
    od;
    return false;
end;


\$SAGE.PlausibleTabCompletionsForSage := function(o)
    local   ops,  opnames;
    ops := Filtered(\$SAGE.OperationsAdmittingFirstArgument(o), op ->
                   \$SAGE.HasAtLeastOneMethodAsFirstArgument(op,o));
    opnames := List(ops, op -> \$SAGE.CleanOperationName(NameFunction(op)));
    return Concatenation(opnames, GLOBAL_FUNCTION_NAMES);
end;

# The log below is for debuging only.
# CAREFUL -- do *not* activate this unless you know
# what you are doing.  E.g., if active and the user doesn't
# have write permission to /tmp (e.g., on OS X),
# then gap will completely fail to work for them.   -- WAS
#
# LogTo("/tmp/gapsage.log");
#
