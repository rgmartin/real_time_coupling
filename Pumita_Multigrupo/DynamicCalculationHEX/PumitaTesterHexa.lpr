program PumitaTesterHexa;

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX}{$IFDEF UseCThreads}
  cthreads,
  {$ENDIF}{$ENDIF}
  Classes, SysUtils, CustApp   ,math,
  reactor_hexa,CAREMHEXA;

const
  NG = 5;
  NTROZOS = 14;
  NCANALES = 61;

type

  { PumitaTesterHexa}

  PumitaTesterHexaProgram = class(TCustomApplication)
  protected
    procedure DoRun; override;
  public
  end;

{ PumitaTesterHexaProgram }

GroupType  = array[1..NG] of real;
   CoreMatrix = array[1..NCANALES,1..NTROZOS] of real;
   CoreFluxMatrix = array[1..NCANALES,1..NTROZOS] of GroupType;
   ChannelType = array[1..NCANALES] of real;
   MemType = array[0..2500000] of byte;
   TString = array [0..256] of char;

var
  TT0Init: real;
  TFinal: real;
  TTime: real;
  DeltaT: real;
  TCPU_no_estacionario: real;
  TCPU_por_ciclo: real;
  TCPU_max_por_ciclo: real;
  TCPU_total: real;

(*******************************************************************************)
 function hacer_strr(a: single; n,d: integer) : string;
 var
   s: string;

 begin
   str (a:n:d,s);
   hacer_strr := s;
   end;   (*  hacer_strr   *)

(*******************************************************************************)
 procedure PosicionarBarras (t: real);    //se simula un scram a partir de t=2
 var
   B: real;
   B01,B02,B09,B11,B13,B03,B08,B10,B12,B07: real;

 begin
   
   if t<10 then
      B:=40.0
   else
       B:=140.0;

   B01 := 70.0;
   B02 := 120.0;
   B09 := B;
   B11 := 0.0;
   B13 := 0.0;
   B03 := 0.0;
   B08 := 0.0;
   B10 := 0.0 ;
   B12 := 0.0;
   B07 := 0.0;

   InsercionBancos  (B01,B02,B09,B11,B13,B03,B08,B10,B12,B07);

   end;

(*******************************************************************************)
 procedure CalculoCinetico;
 var
   NIT: integer;
   TTInit: real;
   TTInitCiclo: real;
   DuracionCiclo: real;
   DuracionMax: real;
   E5,E4,D6,F5,D3,C8,H4,H2,B5,E8,F2,B7,G6,F4,D5,E6,E2,B8,H5,H3,C4,D8,G2,B6,F7,D4,D7,G4:real;
   fBarras: TextFile;

 begin
   TTInit := Time;
   NIT := 0;
   DuracionMax := 0.0;



   while TTime < 0.999999*TFinal do begin
     TTime := TTime + DeltaT;

     TTInitCiclo := Time;

     PosicionarBarras (TTime);

     UnPaso (TTime,DeltaT);
     DuracionCiclo := Time - TTInitCiclo;
     DuracionMax := max(DuracionCiclo,DuracionMax);

     imprimir_salida(TTime,DeltaT);


     inc(NIT);
     end;


   imprimir_salida(TTime,DeltaT);


   //TCPU_no_estacionario := (Time-TTInit)*86400.0;
   //TCPU_por_ciclo := (Time-TTInit)*86400.0/NIT;
   //TCPU_max_por_ciclo := DuracionMax*86400.0;

   end;

(*******************************************************************************)
 procedure EstadoEstacionario;
 var
   n, canal, trozo: integer;
   B01,B02,B09,B11,B13,B03,B08,B10,B12,B07: real;
   cero: real;
   PPMBoro: CoreMatrix;


 begin
   B01 := 70.0;
   B02 := 120.0;
   B09 := 40.0;
   B11 := 0.0;
   B13 := 0.0;
   B03 := 0.0;
   B08 := 0.0;
   B10 := 0.0;
   B12 := 0.0;
   B07 := 0.0;
   cero := 0.0;

   for canal:=1 to NCANALES do
   for trozo:=1 to NTROZOS do
     PPMBoro[canal,trozo] := 0.0;
   Pasar_BORO(PPMBoro);


   InsercionBancos  (B01,B02,B09,B11,B13,B03,B08,B10,B12,B07);
   

   for n:=1 to 2 do begin
     UnPaso(cero,cero);
     end;

   imprimir_salida(cero,cero);
   end;






procedure PumitaTesterHexaProgram.DoRun;
var
  ErrorMsg: String;
begin

  TTime := 0.0;
    TT0Init := Time;
    DeltaT := 5; //5 minutos
    TCPU_no_estacionario := 0.0;
    TCPU_por_ciclo := 0.0;
    TCPU_max_por_ciclo := 0.0;
    TFinal := 1000;//10 horas

    EstadoEstacionario;
    CalculoCinetico;

    TCPU_total := (Time-TT0Init)*86400.0;

    Imprimir_Tiempos_CPU(TCPU_no_estacionario, TCPU_por_ciclo,
                     TCPU_max_por_ciclo,TCPU_total);

    writeln ('Tiempo Total: ' +  hacer_strr (TCPU_total,1,3)) ;
    readln(TCPU_total);
  Terminate;
end;

var
  Application: PumitaTesterHexaProgram;
begin
  Application:=PumitaTesterHexaProgram.Create(nil);
  Application.Title:='PumitaTesterHexaProgram';
  Application.Run;
  Application.Free;
end.

