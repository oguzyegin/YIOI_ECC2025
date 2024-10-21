function [funcName,warning]=funcFind(enteredFunc)
try 
    myFunc = str2func(enteredFunc);
    myFunc(1);
    funcName = enteredFunc;
    warning = 1;
catch
    warning = 1;
    try
        switch enteredFunc
            case 'ln'
                funcName='log';
            case ['log',num2str(enteredFunc(4:end))]
                funcName=['log',num2str(enteredFunc(4:end))];
                warning = 2;
            otherwise
                warning = 0;
                funcName = 'Defined function could not be found in database!'
        end
    catch
        warning = 0;
        funcName = 'Defined function could not be found in database!'
    end
end