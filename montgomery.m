// load "/home/ozbayelif/Development/FIWE/ecm/montgomery.m";
/****************************************************************************/

addM:=function(x1,y1,x2,y2,a,b)
    return b*(y2-y1)^2/(x2-x1)^2-a-x1-x2,(2*x1+x2+a)*(y2-y1)/(x2-x1)-b*(y2-y1)^3/(x2-x1)^3-y1;
end function;

dblM:=function(x1,y1,a,b)
    return b*(3*x1^2+2*a*x1+1)^2/(2*b*y1)^2-a-x1-x1,(2*x1+x1+a)*(3*x1^2+2*a*x1+1)/(2*b*y1)-b*(3*x1^2+2*a*x1+1)^3/(2*b*y1)^3-y1;
end function;

ADDM:=function(X1,Z1,X2,Z2,Xd,Zd)
    A:=X2+Z2;
    B:=X2-Z2;
    C:=X1+Z1;
    D:=X1-Z1;
    DA:=D*A;
    CB:=C*B;
    return Zd*(DA+CB)^2,Xd*(DA-CB)^2;
end function;

DBLM:=function(X1,Z1,A24)
    A:=X1+Z1;
    AA:=A^2;
    B:=X1-Z1;
    BB:=B^2;
    C:=AA-BB;
    return AA*BB,C*(BB+A24*C);
end function;

LADDM:=function(X1,Z1,k,A24)
    R0X:=X1;
    R0Z:=Z1;
    R1X,R1Z:=DBLM(X1,Z1,A24);
    kbit:=Intseq(k,2);
    for i:=#Intseq(k,2)-1 to 1 by -1 do
        if kbit[i] eq 0 then
            R1X,R1Z:=ADDM(R0X,R0Z,R1X,R1Z,X1,Z1);
            R0X,R0Z:=DBLM(R0X,R0Z,A24);
        else
            R0X,R0Z:=ADDM(R0X,R0Z,R1X,R1Z,X1,Z1);
            R1X,R1Z:=DBLM(R1X,R1Z,A24);
        end if;
    end for;
    return R0X,R0Z;
end function;

/****************************************************************************/