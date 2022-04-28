(*$R-*)
 unit interpol;
 interface
   uses sysutils;
 const
   NumXS = 40;
   NG = 5;
 type
   TXS = array[1..NG*NG + 3*NG] of real;

 function Interpolate (P1,P2,P3,P4,P5,P6,P7: real; NTable: integer): TXS;
 function Interpolate_Linear (P1,P2,P3,P4,P5,P6,P7: real; NTable: integer): TXS;
 procedure SetFileName (s: string; NTable,NParm: integer);

 implementation
   uses math, Dialogs;

 const
   SX00: TXS = (0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
   MaxNumTabs = 18;

 type
   TTable2 = array of array of TXS;
   TTable3 = array of TTable2;
   TTable4 = array of TTable3;
   TTable5 = array of TTable4;
   TTable6 = array of TTable5;
   TTable7 = array of TTable6;
   TValues = array of real;

   TAllVars = record
     Parm1Values: TValues;
     Parm2Values: TValues;
     Parm3Values: TValues;
     Parm4Values: TValues;
     Parm5Values: TValues;
     Parm6Values: TValues;
     Parm7Values: TValues;
     table3: TTable3;
     table4: TTable4;
     table5: TTable5;
     table6: TTable6;
     table7: TTable7;
     FileName: string;
     ND1,ND2,ND3,ND4,ND5,ND6,ND7: integer;
     NumDims: integer;
     initialized: boolean;
     end;

 var
   AllVars: array[1..MaxNumTabs] of TAllVars;
   ActualTable: integer;

(*****************************************************************************)
 procedure allocation3;
 var
   i,j,k: integer;

 begin
 with AllVars[ActualTable] do begin
   SetLength(Table3,ND1);

   for i:=0 to pred(ND1) do
     SetLength(Table3[i],ND2);

   for i:=0 to pred(ND1) do
   for j:=0 to pred(ND2) do
     SetLength(Table3[i,j],ND3);

   for i:=0 to pred(ND1) do
   for j:=0 to pred(ND2) do
   for k:=0 to pred(ND3) do
     Table3[i,j,k] := SX00;

   SetLength(Parm1Values,ND1);
   SetLength(Parm2Values,ND2);
   SetLength(Parm3Values,ND3);

   initialized := true;
   end; end;

(*****************************************************************************)
 procedure allocation4;
 var
   i,j,k,m: integer;

 begin
 with AllVars[ActualTable] do begin
   SetLength(Table4,ND1);

   for i:=0 to pred(ND1) do
     SetLength(Table4[i],ND2);

   for i:=0 to pred(ND1) do
   for j:=0 to pred(ND2) do
     SetLength(Table4[i,j],ND3);

   for i:=0 to pred(ND1) do
   for j:=0 to pred(ND2) do
   for k:=0 to pred(ND3) do
     SetLength(Table4[i,j,k],ND4);

   for i:=0 to pred(ND1) do
   for j:=0 to pred(ND2) do
   for k:=0 to pred(ND3) do
   for m:=0 to pred(ND4) do
     Table4[i,j,k,m] := SX00;

   SetLength(Parm1Values,ND1);
   SetLength(Parm2Values,ND2);
   SetLength(Parm3Values,ND3);
   SetLength(Parm4Values,ND4);

   initialized := true;
   end; end;

(*****************************************************************************)
 procedure allocation5;
 var
   i,j,k,m,n: integer;

 begin
 with AllVars[ActualTable] do begin
   SetLength(Table5,ND1);

   for i:=0 to pred(ND1) do
     SetLength(Table5[i],ND2);

   for i:=0 to pred(ND1) do
   for j:=0 to pred(ND2) do
     SetLength(Table5[i,j],ND3);

   for i:=0 to pred(ND1) do
   for j:=0 to pred(ND2) do
   for k:=0 to pred(ND3) do
     SetLength(Table5[i,j,k],ND4);

   for i:=0 to pred(ND1) do
   for j:=0 to pred(ND2) do
   for k:=0 to pred(ND3) do
   for m:=0 to pred(ND4) do
     SetLength(Table5[i,j,k,m],ND5);

   for i:=0 to pred(ND1) do
   for j:=0 to pred(ND2) do
   for k:=0 to pred(ND3) do
   for m:=0 to pred(ND4) do
   for n:=0 to pred(ND5) do
     Table5[i,j,k,m,n] := SX00;

   SetLength(Parm1Values,ND1);
   SetLength(Parm2Values,ND2);
   SetLength(Parm3Values,ND3);
   SetLength(Parm4Values,ND4);
   SetLength(Parm5Values,ND5);

   initialized := true;
   end; end;

(*****************************************************************************)
 procedure allocation6;
 var
   i,j,k,m,n,r: integer;

 begin
 with AllVars[ActualTable] do begin
   SetLength(Table6,ND1);

   for i:=0 to pred(ND1) do
     SetLength(Table6[i],ND2);

   for i:=0 to pred(ND1) do
   for j:=0 to pred(ND2) do
     SetLength(Table6[i,j],ND3);

   for i:=0 to pred(ND1) do
   for j:=0 to pred(ND2) do
   for k:=0 to pred(ND3) do
     SetLength(Table6[i,j,k],ND4);

   for i:=0 to pred(ND1) do
   for j:=0 to pred(ND2) do
   for k:=0 to pred(ND3) do
   for m:=0 to pred(ND4) do
     SetLength(Table6[i,j,k,m],ND5);

   for i:=0 to pred(ND1) do
   for j:=0 to pred(ND2) do
   for k:=0 to pred(ND3) do
   for m:=0 to pred(ND4) do
   for n:=0 to pred(ND5) do
     SetLength(Table6[i,j,k,m],ND6);


   for i:=0 to pred(ND1) do
   for j:=0 to pred(ND2) do
   for k:=0 to pred(ND3) do
   for m:=0 to pred(ND4) do
   for n:=0 to pred(ND5) do
   for r:=0 to pred(ND6) do
     Table6[i,j,k,m,n,r] := SX00;

   SetLength(Parm1Values,ND1);
   SetLength(Parm2Values,ND2);
   SetLength(Parm3Values,ND3);
   SetLength(Parm4Values,ND4);
   SetLength(Parm5Values,ND5);
   SetLength(Parm6Values,ND6);

   initialized := true;
   end; end;

(*****************************************************************************)
 procedure allocation7;
 var
   i,j,k,m,n,r,s: integer;

 begin
 with AllVars[ActualTable] do begin
   SetLength(Table7,ND1);

   for i:=0 to pred(ND1) do
     SetLength(Table7[i],ND2);

   for i:=0 to pred(ND1) do
   for j:=0 to pred(ND2) do
     SetLength(Table7[i,j],ND3);

   for i:=0 to pred(ND1) do
   for j:=0 to pred(ND2) do
   for k:=0 to pred(ND3) do
     SetLength(Table7[i,j,k],ND4);

   for i:=0 to pred(ND1) do
   for j:=0 to pred(ND2) do
   for k:=0 to pred(ND3) do
   for m:=0 to pred(ND4) do
     SetLength(Table7[i,j,k,m],ND5);

   for i:=0 to pred(ND1) do
   for j:=0 to pred(ND2) do
   for k:=0 to pred(ND3) do
   for m:=0 to pred(ND4) do
   for n:=0 to pred(ND5) do
     SetLength(Table7[i,j,k,m,n],ND6);

   for i:=0 to pred(ND1) do
   for j:=0 to pred(ND2) do
   for k:=0 to pred(ND3) do
   for m:=0 to pred(ND4) do
   for n:=0 to pred(ND5) do
   for r:=0 to pred(ND6) do
     SetLength(Table7[i,j,k,m,n,r],ND7);

   for i:=0 to pred(ND1) do
   for j:=0 to pred(ND2) do
   for k:=0 to pred(ND3) do
   for m:=0 to pred(ND4) do
   for n:=0 to pred(ND5) do
   for r:=0 to pred(ND6) do
   for s:=0 to pred(ND7) do
     Table7[i,j,k,m,n,r,s] := SX00;

   SetLength(Parm1Values,ND1);
   SetLength(Parm2Values,ND2);
   SetLength(Parm3Values,ND3);
   SetLength(Parm4Values,ND4);
   SetLength(Parm5Values,ND5);
   SetLength(Parm6Values,ND6);
   SetLength(Parm7Values,ND7);

   initialized := true;
   end; end;

 (*****************************************************************************)
 procedure SetFileName (s: string; NTable,NParm: integer);
 begin
 ActualTable := NTable;

 with AllVars[ActualTable] do begin
   FileName := s;

   if not FileExists(s) then begin
     ShowMessage('No existe el archivo "' + s + '"');
     halt(0);
     end;

   NumDims := NParm+1;
   end; end;

 (*****************************************************************************)
 procedure ReadTable;
 var
   ii: integer;
   i,j,k,m,n,q,s: integer;
   F: text;
   R: TXS;
   RR: array[1..100] of real;
   NDX: integer;

 begin
 with AllVars[ActualTable] do begin
   AssignFile(F,AllVars[ActualTable].FileName);
   Reset(F);

   for ii:=1 to 100 do RR[ii] := 0.0;

   for ii:=1 to 10 do begin Read(F,RR[ii]);  end;

   ND1 := trunc(RR[1]);
   ND2 := trunc(RR[2]);
   ND3 := trunc(RR[3]);
   ND4 := trunc(RR[4]);
   ND5 := trunc(RR[5]);
   ND6 := trunc(RR[6]);
   ND7 := trunc(RR[7]);

   NDX  := 7 - ord(ND7 = 0) - ord(ND6 = 0) - ord(ND5 = 0) - ord(ND4 = 0);

   if NDX <> NumDims then begin
     ShowMessage ('El número de dimensiones de la tabla leída no coincide '
                      + 'con el declarado');
     Halt(0);
     end;

   case NumDims of
     3: Allocation3;
     4: Allocation4;
     5: Allocation5;
     6: Allocation6;
     7: Allocation7;
     end;

   readln(F);
   for ii:=1 to ND1 do begin Read(F,RR[ii]); end;
   for i:=1 to ND1 do Parm1Values[i-1] := RR[i];

   readln(F);
   for ii:=1 to ND2 do Read(F,RR[ii]);
   for i:=1 to ND2 do Parm2Values[i-1] := RR[i];

   readln(F);
   for ii:=1 to ND3 do Read(F,RR[ii]);
   for i:=1 to ND3 do Parm3Values[i-1] := RR[i];

   if NumDims > 3 then begin
     readln(F);
     for ii:=1 to ND4 do Read(F,RR[ii]);
     for i:=1 to ND4 do Parm4Values[i-1] := RR[i];
     end;

   if NumDims > 4 then begin
     readln(F);
     for ii:=1 to ND5 do Read(F,RR[ii]);
     for i:=1 to ND5 do Parm5Values[i-1] := RR[i];
     end;

   if NumDims > 5 then begin
     readln(F);
     for ii:=1 to ND6 do Read(F,RR[ii]);
     for i:=1 to ND6 do Parm6Values[i-1] := RR[i];
     end;

   if NumDims > 6 then begin
     readln(F);
     for ii:=1 to ND7 do Read(F,RR[ii]);
     for i:=1 to ND7 do Parm7Values[i-1] := RR[i];
     end;


   case Numdims of
     3: for i:=0 to pred(ND1) do
        for j:=0 to pred(ND2) do
        for k:=0 to pred(ND3) do begin
          readln(F);
          for ii:=1 to NumXS do Read(F,R[ii]);
          Table3[i,j,k] := R;
          end;

     4: for i:=0 to pred(ND1) do
        for j:=0 to pred(ND2) do
        for k:=0 to pred(ND3) do
        for m:=0 to pred(ND4) do begin
          readln(F);
          for ii:=1 to NumXS do Read(F,R[ii]);
          Table4[i,j,k,m] := R;
          end;

     5: for i:=0 to pred(ND1) do
        for j:=0 to pred(ND2) do
        for k:=0 to pred(ND3) do
        for m:=0 to pred(ND4) do
        for n:=0 to pred(ND5) do begin
          readln(F);
          for ii:=1 to NumXS do Read(F,R[ii]);
          Table5[i,j,k,m,n] := R;
          end;

     6: for i:=0 to pred(ND1) do
        for j:=0 to pred(ND2) do
        for k:=0 to pred(ND3) do
        for m:=0 to pred(ND4) do
        for n:=0 to pred(ND5) do
        for q:=0 to pred(ND6) do begin
          readln(F);
          for ii:=1 to NumXS do Read(F,R[ii]);
          Table6[i,j,k,m,n,q] := R;
          end;

     7: for i:=0 to pred(ND1) do
        for j:=0 to pred(ND2) do
        for k:=0 to pred(ND3) do
        for m:=0 to pred(ND4) do
        for n:=0 to pred(ND5) do
        for q:=0 to pred(ND6) do
        for s:=0 to pred(ND7) do begin
          readln(F);
          for ii:=1 to NumXS do Read(F,R[ii]);
          Table7[i,j,k,m,n,q,s] := R;
          end;

     end;

   Close(F);
   end; end;

 (*****************************************************************************)
 procedure ReadTableBinary;
 var
   i,j,k,m,n,q,s: integer;
   F: file;
   R: TXS;
   RR: array[1..100] of real;
   NDX: integer;

 begin
 with AllVars[ActualTable] do begin
   AssignFile(F,FileName);
   Reset(F,1);

   BlockRead(F,R,SizeOf(R));
   ND1 := trunc(R[1]);
   ND2 := trunc(R[2]);
   ND3 := trunc(R[3]);
   ND4 := trunc(R[4]);
   ND5 := trunc(R[5]);
   ND6 := trunc(R[6]);
   ND7 := trunc(R[7]);

   NDX  := 7 - ord(ND7 = 0) - ord(ND6 = 0) - ord(ND5 = 0) - ord(ND4 = 0);

   if NDX <> NumDims then begin
     ShowMessage ('El número de dimensiones de la tabla leída no coincide '
                      + 'con el declarado');
     Halt(0);
     end;

   case NumDims of
     3: Allocation3;
     4: Allocation4;
     5: Allocation5;
     6: Allocation6;
     7: Allocation7;
     end;

   BlockRead(F,RR,ND1*SizeOf(real));
   for i:=1 to ND1 do Parm1Values[i-1] := RR[i];

   BlockRead(F,RR,ND2*SizeOf(real));
   for i:=1 to ND2 do Parm2Values[i-1] := RR[i];

   BlockRead(F,RR,ND3*SizeOf(real));
   for i:=1 to ND3 do Parm3Values[i-1] := RR[i];

   if NumDims > 3 then begin
     BlockRead(F,RR,ND4*SizeOf(real));
     for i:=1 to ND4 do Parm4Values[i-1] := RR[i];
     end;

   if NumDims > 4 then begin
     BlockRead(F,RR,ND5*SizeOf(real));
     for i:=1 to ND5 do Parm5Values[i-1] := RR[i];
     end;

   if NumDims > 5 then begin
     BlockRead(F,RR,ND6*SizeOf(real));
     for i:=1 to ND6 do Parm6Values[i-1] := RR[i];
     end;

   if NumDims > 6 then begin
     BlockRead(F,RR,ND7*SizeOf(real));
     for i:=1 to ND7 do Parm7Values[i-1] := RR[i];
     end;


   case Numdims of
     3: for i:=0 to pred(ND1) do
        for j:=0 to pred(ND2) do
        for k:=0 to pred(ND3) do begin
          BlockRead(F,R,SizeOf(R));
          Table3[i,j,k] := R;
          end;

     4: for i:=0 to pred(ND1) do
        for j:=0 to pred(ND2) do
        for k:=0 to pred(ND3) do
        for m:=0 to pred(ND4) do begin
          BlockRead(F,R,SizeOf(R));
          Table4[i,j,k,m] := R;
          end;

     5: for i:=0 to pred(ND1) do
        for j:=0 to pred(ND2) do
        for k:=0 to pred(ND3) do
        for m:=0 to pred(ND4) do
        for n:=0 to pred(ND5) do begin
          BlockRead(F,R,SizeOf(R));
          Table5[i,j,k,m,n] := R;
          end;

     6: for i:=0 to pred(ND1) do
        for j:=0 to pred(ND2) do
        for k:=0 to pred(ND3) do
        for m:=0 to pred(ND4) do
        for n:=0 to pred(ND5) do
        for q:=0 to pred(ND6) do begin
          BlockRead(F,R,SizeOf(R));
          Table6[i,j,k,m,n,q] := R;
          end;

     7: for i:=0 to pred(ND1) do
        for j:=0 to pred(ND2) do
        for k:=0 to pred(ND3) do
        for m:=0 to pred(ND4) do
        for n:=0 to pred(ND5) do
        for q:=0 to pred(ND6) do
        for s:=0 to pred(ND7) do begin
          BlockRead(F,R,SizeOf(R));
          Table7[i,j,k,m,n,q,s] := R;
          end;

     end;

   Close(F);
   end; end;

(*****************************************************************************)
 procedure locateQ (T: array of real; val: real;
                   var i0,i1,i2: integer; var C0,C1,C2: real);
 var
   k, N, i: integer;
   x1,x2,x3: real;

 begin
   N := High(T);

   k := 1;
   while (k <= N) and ((val-T[k]) >= 0) do inc(k);
//   DO k:=1 TO N WHILE ((val-T(k,1)) >= 0); END;
     k := k - 1;
     i := min(max(k,1),N-1);

     i0 := i-1;
     i1 := i;
     i2 := i+1;

     x1 := T[i0];
     x2 := T[i1];
     x3 := T[i2];

     C0 := (val-x2)*(val-x3)/((x1-x2)*(x1-x3));
     C1 := (val-x1)*(val-x3)/((x2-x1)*(x2-x3));
     C2 := (val-x1)*(val-x2)/((x3-x1)*(x3-x2));
     end;

(*****************************************************************************)
 procedure locate (T: array of real; val: real; var i1,i2: integer; var C1,C2: real);
 var
   i, N: integer;

 begin
   N := High(T);

   if N = 0 then begin
     i1 := 0;
     i2 := 0;
     C1 := 1.0;
     C2 := 0.0;
     end

   else if N = 1 then begin
     i1 := 0;
     i2 := 1;
     C2 := (val - T[0])/(T[1] - T[0]);
     C1 := 1.0 - C2;
     end

   else begin
     i := 0;

     while (T[i] < val) and (i <= N) do inc(i);

     i := Max(1,min(i,N));
     i1 := i-1;
     i2 := i;
     C2 := (val - T[i1])/(T[i2] - T[i1]);
     C1 := 1.0 - C2;
     end;

   end;

(*****************************************************************************)
 function Evaluate2 (var Table2: TTable2; k1,m1,k2,m2: integer;
                     ck1,cm1,ck2,cm2: real): TXS;
 var
   n: integer;

 begin

   for n:=1 to NumXS do
   Result[n] :=    ck1*(cm1*Table2[k1,m1,n] +
                        cm2*Table2[k1,m2,n]) +
                   ck2*(cm1*Table2[k2,m1,n] +
                        cm2*Table2[k2,m2,n]);

   end;

(*****************************************************************************)
 function Evaluate3 (var Table3: TTable3; j1,k1,m1,j2,k2,m2: integer;
                     cj1,ck1,cm1,cj2,ck2,cm2: real): TXS;
 var
   n: integer;

 begin

   for n:=1 to NumXS do
   Result[n] :=  (cj1*(ck1*(cm1*Table3[j1,k1,m1,n] +
                       cm2*Table3[j1,k1,m2,n]) +
                       ck2*(cm1*Table3[j1,k2,m1,n] +
                       cm2*Table3[j1,k2,m2,n])) +
                  cj2*(ck1*(cm1*Table3[j2,k1,m1,n] +
                       cm2*Table3[j2,k1,m2,n]) +
                       ck2*(cm1*Table3[j2,k2,m1,n] +
                       cm2*Table3[j2,k2,m2,n])));

   end;

(*****************************************************************************)
 function Evaluate4 (var Table4: TTable4; i1,j1,k1,m1,i2,j2,k2,m2: integer;
                     ci1,cj1,ck1,cm1,ci2,cj2,ck2,cm2: real): TXS;
 var
   n: integer;

 begin

   for n:=1 to NumXS do
   Result[n] :=
            ci1*(cj1*(ck1*(cm1*Table4[i1,j1,k1,m1,n] +
                           cm2*Table4[i1,j1,k1,m2,n]) +
                      ck2*(cm1*Table4[i1,j1,k2,m1,n] +
                           cm2*Table4[i1,j1,k2,m2,n])) +
                 cj2*(ck1*(cm1*Table4[i1,j2,k1,m1,n] +
                           cm2*Table4[i1,j2,k1,m2,n]) +
                      ck2*(cm1*Table4[i1,j2,k2,m1,n] +
                           cm2*Table4[i1,j2,k2,m2,n]))) +

            ci2*(cj1*(ck1*(cm1*Table4[i2,j1,k1,m1,n] +
                           cm2*Table4[i2,j1,k1,m2,n]) +
                      ck2*(cm1*Table4[i2,j1,k2,m1,n] +
                           cm2*Table4[i2,j1,k2,m2,n])) +
                 cj2*(ck1*(cm1*Table4[i2,j2,k1,m1,n] +
                           cm2*Table4[i2,j2,k1,m2,n]) +
                      ck2*(cm1*Table4[i2,j2,k2,m1,n] +
                           cm2*Table4[i2,j2,k2,m2,n])));

   end;

(*****************************************************************************)
 function Evaluate5 (var Table5: TTable5; i1,j1,k1,m1,n1,i2,j2,k2,m2,n2: integer;
                     ci1,cj1,ck1,cm1,cn1,ci2,cj2,ck2,cm2,cn2: real): TXS;
 var
   n: integer;
   RR1,RR2: TXS;

 begin
    RR1 := Evaluate4 (Table5[i1],j1,k1,m1,n1,j2,k2,m2,n2,
                               cj1,ck1,cm1,cn1,cj2,ck2,cm2,cn2);
    RR2 := Evaluate4 (Table5[i2],j1,k1,m1,n1,j2,k2,m2,n2,
                               cj1,ck1,cm1,cn1,cj2,ck2,cm2,cn2);

    for n:=1 to NumXS do
      Result[n] := ci1*RR1[n] + ci2*RR2[n];

   end;

(*****************************************************************************)
 function Evaluate6 (var Table6: TTable6; i1,j1,k1,m1,n1,r1,i2,j2,k2,m2,n2,r2: integer;
                     ci1,cj1,ck1,cm1,cn1,cr1,ci2,cj2,ck2,cm2,cn2,cr2: real): TXS;
 var
   n: integer;
   RR1,RR2: TXS;

 begin
    RR1 := Evaluate5 (Table6[i1],j1,k1,m1,n1,r1,j2,k2,m2,n2,r2,
                               cj1,ck1,cm1,cn1,cr1,cj2,ck2,cm2,cn2,cr2);
    RR2 := Evaluate5 (Table6[i2],j1,k1,m1,n1,r1,j2,k2,m2,n2,r2,
                               cj1,ck1,cm1,cn1,cr1,cj2,ck2,cm2,cn2,cr2);

    for n:=1 to NumXS do
      Result[n] := ci1*RR1[n] + ci2*RR2[n];

   end;

(*****************************************************************************)
 function Interpolate_Linear (P1,P2,P3,P4,P5,P6,P7: real; NTable: integer): TXS;
 var
   n: integer;
//   i0: integer;
   i1,j1,k1,m1,n1,r1,s1: integer;
   i2,j2,k2,m2,n2,r2,s2: integer;
//   ci0: real;
   ci1,cj1,ck1,cm1,cn1,cr1,cs1: real;
   ci2,cj2,ck2,cm2,cn2,cr2,cs2: real;
//   RR0,
   RR1,RR2,R: TXS;

 begin
 ActualTable := NTable;

 with AllVars[ActualTable] do begin
   if not initialized then
     if pos('.DAT',FileName) + pos('.dat',FileName) = 0 then ReadTable
                                                        else ReadTableBinary;

//   ci0 := 0.0; i0 := 0; for n:=1 to NumXS do RR0[n] := 0.0;
   locate (Parm1Values,P1,i1,i2,ci1,ci2);
   locate (Parm2Values,P2,j1,j2,cj1,cj2);
   locate (Parm3Values,P3,k1,k2,ck1,ck2);

   if NumDims > 3 then
     locate (Parm4Values,P4,m1,m2,cm1,cm2);

   if NumDims > 4 then
     locate (Parm5Values,P5,n1,n2,cn1,cn2);

   if NumDims > 5 then
     locate (Parm6Values,P6,r1,r2,cr1,cr2);

   if NumDims > 6 then
     locate (Parm7Values,P7,s1,s2,cs1,cs2);

   case NumDims of
     3: begin
//          RR0 := evaluate2(Table3[i0],j1,k1,j2,k2,cj1,ck1,cj2,ck2);
          RR1 := evaluate2(Table3[i1],j1,k1,j2,k2,cj1,ck1,cj2,ck2);
          RR2 := evaluate2(Table3[i2],j1,k1,j2,k2,cj1,ck1,cj2,ck2);
          end;

     4: begin
//          RR0 := evaluate3(Table4[i0],j1,k1,m1,j2,k2,m2,cj1,ck1,cm1,cj2,ck2,cm2);
          RR1 := evaluate3(Table4[i1],j1,k1,m1,j2,k2,m2,cj1,ck1,cm1,cj2,ck2,cm2);
          RR2 := evaluate3(Table4[i2],j1,k1,m1,j2,k2,m2,cj1,ck1,cm1,cj2,ck2,cm2);
          end;

     5: begin
//          RR0 := evaluate4(Table5[i0],j1,k1,m1,n1,j2,k2,m2,n2,
//                                   cj1,ck1,cm1,cn1,cj2,ck2,cm2,cn2);
          RR1 := evaluate4(Table5[i1],j1,k1,m1,n1,j2,k2,m2,n2,
                                   cj1,ck1,cm1,cn1,cj2,ck2,cm2,cn2);
          RR2 := evaluate4(Table5[i2],j1,k1,m1,n1,j2,k2,m2,n2,
                                   cj1,ck1,cm1,cn1,cj2,ck2,cm2,cn2);
          end;

     6: begin
//          RR0 := evaluate5(Table6[i0],j1,k1,m1,n1,r1,j2,k2,m2,n2,r2,
//                                   cj1,ck1,cm1,cn1,cr1,cj2,ck2,cm2,cn2,cr2);
          RR1 := evaluate5(Table6[i1],j1,k1,m1,n1,r1,j2,k2,m2,n2,r2,
                                   cj1,ck1,cm1,cn1,cr1,cj2,ck2,cm2,cn2,cr2);
          RR2 := evaluate5(Table6[i2],j1,k1,m1,n1,r1,j2,k2,m2,n2,r2,
                                   cj1,ck1,cm1,cn1,cr1,cj2,ck2,cm2,cn2,cr2);
          end;

     7: begin
//          RR0 := evaluate6(Table7[i0],j1,k1,m1,n1,r1,s1,j2,k2,m2,n2,r2,s2,
//                             cj1,ck1,cm1,cn1,cr1,cs1,cj2,ck2,cm2,cn2,cr2,cs2);
          RR1 := evaluate6(Table7[i1],j1,k1,m1,n1,r1,s1,j2,k2,m2,n2,r2,s2,
                             cj1,ck1,cm1,cn1,cr1,cs1,cj2,ck2,cm2,cn2,cr2,cs2);
          RR2 := evaluate6(Table7[i2],j1,k1,m1,n1,r1,s1,j2,k2,m2,n2,r2,s2,
                             cj1,ck1,cm1,cn1,cr1,cs1,cj2,ck2,cm2,cn2,cr2,cs2);
          end;

      end;


   for n:=1 to NumXS do
//     R[n] := ci0*RR0[n] + ci1*RR1[n] + ci2*RR2[n];
     R[n] := ci1*RR1[n] + ci2*RR2[n];

   for n:=1 to NG do
     R[n] := 1.0/(3.0*R[n]);

   Result := R;
   end; end;

(*****************************************************************************)
 function Interpolate (P1,P2,P3,P4,P5,P6,P7: real; NTable: integer): TXS;
 var
   n: integer;
   i0: integer;
   i1,j1,k1,m1,n1,r1,s1: integer;
   i2,j2,k2,m2,n2,r2,s2: integer;
   ci0: real;
   ci1,cj1,ck1,cm1,cn1,cr1,cs1: real;
   ci2,cj2,ck2,cm2,cn2,cr2,cs2: real;
   RR0,RR1,RR2,R: TXS;

 begin
 ActualTable := NTable;

 with AllVars[ActualTable] do begin
   if not initialized then
     if pos('.DAT',FileName) + pos('.dat',FileName) = 0 then ReadTable
                                                        else ReadTableBinary;

   ci0 := 0.0; i0 := 0; for n:=1 to NumXS do RR0[n] := 0.0;
   locateQ (Parm1Values,P1,i0,i1,i2, ci0,ci1,ci2);
   locate (Parm2Values,P2,j1,j2,cj1,cj2);
   locate (Parm3Values,P3,k1,k2,ck1,ck2);

   if NumDims > 3 then
     locate (Parm4Values,P4,m1,m2,cm1,cm2);

   if NumDims > 4 then
     locate (Parm5Values,P5,n1,n2,cn1,cn2);

   if NumDims > 5 then
     locate (Parm6Values,P6,r1,r2,cr1,cr2);

   if NumDims > 6 then
     locate (Parm7Values,P7,s1,s2,cs1,cs2);

   case NumDims of
     3: begin
          RR0 := evaluate2(Table3[i0],j1,k1,j2,k2,cj1,ck1,cj2,ck2);
          RR1 := evaluate2(Table3[i1],j1,k1,j2,k2,cj1,ck1,cj2,ck2);
          RR2 := evaluate2(Table3[i2],j1,k1,j2,k2,cj1,ck1,cj2,ck2);
          end;

     4: begin
          RR0 := evaluate3(Table4[i0],j1,k1,m1,j2,k2,m2,cj1,ck1,cm1,cj2,ck2,cm2);
          RR1 := evaluate3(Table4[i1],j1,k1,m1,j2,k2,m2,cj1,ck1,cm1,cj2,ck2,cm2);
          RR2 := evaluate3(Table4[i2],j1,k1,m1,j2,k2,m2,cj1,ck1,cm1,cj2,ck2,cm2);
          end;

     5: begin
          RR0 := evaluate4(Table5[i0],j1,k1,m1,n1,j2,k2,m2,n2,
                                   cj1,ck1,cm1,cn1,cj2,ck2,cm2,cn2);
          RR1 := evaluate4(Table5[i1],j1,k1,m1,n1,j2,k2,m2,n2,
                                   cj1,ck1,cm1,cn1,cj2,ck2,cm2,cn2);
          RR2 := evaluate4(Table5[i2],j1,k1,m1,n1,j2,k2,m2,n2,
                                   cj1,ck1,cm1,cn1,cj2,ck2,cm2,cn2);
          end;

     6: begin
          RR0 := evaluate5(Table6[i0],j1,k1,m1,n1,r1,j2,k2,m2,n2,r2,
                                   cj1,ck1,cm1,cn1,cr1,cj2,ck2,cm2,cn2,cr2);
          RR1 := evaluate5(Table6[i1],j1,k1,m1,n1,r1,j2,k2,m2,n2,r2,
                                   cj1,ck1,cm1,cn1,cr1,cj2,ck2,cm2,cn2,cr2);
          RR2 := evaluate5(Table6[i2],j1,k1,m1,n1,r1,j2,k2,m2,n2,r2,
                                   cj1,ck1,cm1,cn1,cr1,cj2,ck2,cm2,cn2,cr2);
          end;

     7: begin
          RR0 := evaluate6(Table7[i0],j1,k1,m1,n1,r1,s1,j2,k2,m2,n2,r2,s2,
                             cj1,ck1,cm1,cn1,cr1,cs1,cj2,ck2,cm2,cn2,cr2,cs2);
          RR1 := evaluate6(Table7[i1],j1,k1,m1,n1,r1,s1,j2,k2,m2,n2,r2,s2,
                             cj1,ck1,cm1,cn1,cr1,cs1,cj2,ck2,cm2,cn2,cr2,cs2);
          RR2 := evaluate6(Table7[i2],j1,k1,m1,n1,r1,s1,j2,k2,m2,n2,r2,s2,
                             cj1,ck1,cm1,cn1,cr1,cs1,cj2,ck2,cm2,cn2,cr2,cs2);
          end;

      end;


   for n:=1 to NumXS do
     R[n] := ci0*RR0[n] + ci1*RR1[n] + ci2*RR2[n];

//   for n:=1 to NG do
//     R[n] := 1.0/(3.0*R[n]);

   Result := R;
   end; end;

(*****************************************************************************)

 begin
   for ActualTable:=1 to MaxNumTabs do
   with AllVars[ActualTable] do begin
     FileName := '';
     initialized := false;
     ND1 := 0; ND2 := 0; ND3 := 0;
     ND4 := 0; ND5 := 0; ND6 := 0; ND7 := 0;
     NumDims := 0;
     end;

   ActualTable := 1;
   end.

