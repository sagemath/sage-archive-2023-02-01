# myextension.py
from IPython.core.inputtransformer import CoroutineInputTransformer
import token
from tokenize import generate_tokens, TokenError
from ast import literal_eval

def string_decorator_defaults(func):
    def my_wrap(*args,**kwds):
        if len(kwds)==0 and len(args)==1 and isinstance(args[0],basestring):
            # call without parentheses
            return func(*args)
        else:
            return lambda f: func(f, *args, **kwds)
    return my_wrap

def match(g, begin):
    #find closing paren; that will be func_end
    nest={'(':0,')':0,'[':0,']':0}
    nest[begin]+=1
    code=token.OP
    while (nest['(']>nest[')'] or nest['[']>nest[']']) and code!=token.ENDMARKER:
        code,value,start,end,_ = g.next()
        if code==token.OP and value in "()[]":
            nest[value]+=1
    if nest['(']!=nest[')'] and nest['[']!=nest[']']:
        # couldn't find a matching closing delimiter
        raise ValueError('Unmatched delimiters')
    return end

def callable(g, code, value, callable_start, callable_end):
    """
    returns the start and end position of the next full callable expression.

    A callable expression consists of a name followed by .name, matched () or matched [].
    """
    if code!=token.NAME:
        raise ValueError("No callable found")
    code,value,start,end,_ = g.next()
    while code==token.OP:
        if value=='.':
            # dot must be followed by a name
            code,value,start,end,_ = g.next()
            if code != token.NAME:
                raise ValueError("No callable found")
            else:
                callable_end = end
        elif value in '([':
            try:
                callable_end = match(g, value)
            except ValueError:
                raise ValueError("No callable found")
        code,value,start,end,_ = g.next()

    return (callable_start, callable_end), (code, value, start, end)



@CoroutineInputTransformer.wrap
def stringdecorator(end_on_blank_line=False):
    """Captures & transforms cell magics.
    
    After a cell magic is started, this stores up any lines it gets until it is
    reset (sent None).
    """
    line = ''
    while True:
        line = (yield line)
        if (not line) or (not line.startswith('%')):
            continue
        g = generate_tokens(iter([line]).next)
        try:
            code,value,start,end,_ = g.next()
            if not (code==token.OP and value == '%'):
                # this check is redundant with the .startswith above
                # but we'll do it anyway just to make sure the %
                # is recognized as an OP
                continue
            code,value,start,end,_ = g.next()
            cell_decorator=False
            end_marker=None
            if code==token.OP and value == '%':
                cell_decorator=True
                # maybe we should just pass %% on to IPython cell decorators
                # and ignore them here?  That's probably too confusing.
                code,value,start,end,_ = g.next()
            # optionally have a string end marker
            if code==token.STRING:
                end_marker = literal_eval(value)
                code,value,start,end,_ = g.next()
            try:
                (func_start, func_end), (code,value,start,end) = callable(g, code,value,start,end)
                func_start=func_start[1]
                func_end=func_end[1]
            except ValueError:
                # ignore this line if we couldn't get a callable
                continue
            decorator = line[func_start:func_end]
            string = line[func_end:].strip()
            if len(string)==0 or cell_decorator:
                body = []
                if len(string)>0: body.append(string)
                line = (yield None)
                while line is not None:
                    if end_marker is not None:
                        if line.strip()==end_marker:
                            break
                    elif end_on_blank_line and line.strip()=='':
                        break
                    body.append(line)
                    line = (yield None)
                string = '\n'.join(body)
            line = "%s(%r)"%(decorator,string)
        except TokenError:
            # something went really wrong
            continue
        print "transformed to %r"%line

def load_ipython_extension(ipython):
    # The `ipython` argument is the currently active `InteractiveShell`
    # instance, which can be used in any way. This allows you to register
    # new magics or aliases, for example.
    ipython.input_splitter.physical_line_transforms.prepend(stringdecorator(end_on_blank_line=True))
    ipython.input_transformer_manager.physical_line_transforms.prepend(stringdecorator(end_on_blank_line=False))

def unload_ipython_extension(ipython):
    # If you want your extension to be unloadable, put that logic here.
    pass

