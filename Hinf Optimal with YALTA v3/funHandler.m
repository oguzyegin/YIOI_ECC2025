function [hNum]=funHandler(hNum,entered)
minus=-1;
if strcmp(entered(1),'-')
    entered = entered(2:end);
    minus=1;
end
if isnan(hNum)
    a = strfind(entered,'(');
    b = strfind(entered,')');
    mul = strfind(entered,'*');
    if mul<a-1
        index=mul+1;
        multip=str2double(entered(1:mul-1));
    else
        index=1;
        multip=1;
    end
    funcName = entered(index:a-1);
    parameter = entered(a+1:b-1);
    [desFunc,bool] = funcFind(funcName);
    if bool==1
        desFunc = str2func(desFunc);
        hNum = multip*minus*desFunc(str2double(parameter));
    elseif bool==2
        desLogBase = str2double(desFunc(4:end));
        desFunc = str2func(desFunc(1:3));
        hNum = multip*minus*desFunc(str2double(parameter))/desFunc(desLogBase);
    else
        if ~isempty(strfind(entered,'/'))
            div_ind = strfind(entered,'/');
            hNum = str2double(entered(1:div_ind-1))/str2double(entered(div_ind+1:end));
            hNum
        else
            msgbox(desFunc,'Error');
        end
    end
end