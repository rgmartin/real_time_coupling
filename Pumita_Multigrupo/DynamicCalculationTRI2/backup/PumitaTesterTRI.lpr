program PumitaTesterTRI;

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX}{$IFDEF UseCThreads}
  cthreads,
  {$ENDIF}{$ENDIF}
  Classes, SysUtils, CustApp,  Dialogs,math,
  reactor,CAREMTRI;

const
  NG = 5;
  NTROZOS = 14;
  NCANALES = 61;


type

  {PumitaTesterTRI }

  pumitaTester = class(TCustomApplication)
  protected
    procedure DoRun; override;
  public
  end;

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
 procedure PosicionarBarras (t: real);    //
 var
   B: real;
   B01,B02,B09,B11,B13,B03,B08,B10,B12,B07: real;

 begin
    //Distribución critica
   B01 := 70.0;    B02 :=120 ; B09 := 55.0; B11 := 0.0;   B13 := 0.0;
   B03 := 0.0;   B08 := 0;   B10 := 0.0 ;   B12 :=00.0;   B07 := 0.0;
   {
   //Reactividad en exceso.. extraer todas las barras
   B01 := 0.0;    B02 :=0 ; B09 := 0; B11 := 0.0;   B13 := 0.0;
   B03 := 0.0;   B08 := 0;   B10 := 0.0 ;   B12 :=00.0;   B07 := 0.0;

   //Reactividad introducida por el SAC
   B01 := 140.0;    B02 :=140.0 ;  B03 := 140.0; B09 := 140.0; B11 := 140.0;   B13 := 140.0;
    B08 := 0;   B10 := 0.0 ;   B12 :=00.0;   B07 := 0.0;

   //Reactividad introducida por el SER
   B01 := 0.0;    B02 :=0.0 ;  B03 := 0.0; B09 := 0.0; B11 := 0.0;   B13 := 0.0;      //SAC
    B08 :=140.0;   B10 :=140.0 ;   B12 :=140.0;                                          //SER
    B07 := 0.0;                                                                       //Resguardo (falta B16)

    }


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

   TCPU_no_estacionario := (Time-TTInit)*86400.0;
   TCPU_por_ciclo := (Time-TTInit)*86400.0/NIT;
   TCPU_max_por_ciclo := DuracionMax*86400.0;

   end;

(*******************************************************************************)
 procedure EstadoEstacionario;
 var
   n, canal, trozo: integer;
   cero:real;
   PPMBoro: CoreMatrix;


 begin

   for canal:=1 to NCANALES do
   for trozo:=1 to NTROZOS do
     PPMBoro[canal,trozo] := 0.0;
   Pasar_BORO(PPMBoro);


   cero:=0.0;
   PosicionarBarras(0);
   for n:=1 to 1 do begin
     UnPaso(cero,cero);
     imprimir_salida(cero,cero);
     end;
   end;

procedure pumitaTester.DoRun;

begin



  //las siguientes lineas pertenecian al miniCAREM original

    TTime := 0;
    TT0Init := Time;
    DeltaT := 1; //1 segnudo
    TCPU_no_estacionario := 0.0;
    TCPU_por_ciclo := 0.0;
    TCPU_max_por_ciclo := 0.0;
    TFinal := 0; //2 minutos
    
    EstadoEstacionario;

    if TFINAL > 0.0 then
       CalculoCinetico;


    TCPU_total := (Time-TT0Init)*86400.0;
    Imprimir_Tiempos_CPU(TCPU_no_estacionario, TCPU_por_ciclo,
                     TCPU_max_por_ciclo,TCPU_total);

    writeln ('Tiempo Total: ' +  hacer_strr (TCPU_total,1,3)) ;
    readln(TCPU_total);
  Terminate;
end;

var
  Application: pumitaTester;
begin
  Application:=pumitaTester.Create(nil);
  Application.Title:='pumitaTester';
  Application.Run;
  Application.Free;
end.

