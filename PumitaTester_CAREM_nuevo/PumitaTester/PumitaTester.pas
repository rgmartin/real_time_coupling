//{$DEFINE UseDLL}
{$DEFINE UseHEXA}
program PumitaTester;
uses math,Dialogs,



{$IFNDEF UseDLL}
  AlReactor,
  {$IFDEF UseHEXA}
   CAREMHEXA,
  {$ELSE}
   CAREMTRI,
  {$ENDIF}
{$ENDIF}
 sysutils;
  
const
  NG = 5;
  NTROZOS = 14;
  NCANALES = 61;

{$IFDEF UseDLL}
{$IFDEF UseHEXA}
  archivo_dll='reactor_hexa.dll' ;
{$ELSE}
  archivo_dll = 'reactor.dll';
{$ENDIF}
  type
    GroupType  = array[1..NG] of real;
    CoreMatrix = array[1..NCANALES,1..NTROZOS] of real;
    CoreFluxMatrix = array[1..NCANALES,1..NTROZOS] of GroupType;
    ChannelType = array[1..NCANALES] of real;
    MemType = array[0..2500000] of byte;
    TString = array [0..256] of char;
  procedure InsercionBancos
     (var B01,B02,B09,B11,B13,B03,B08,B10,B12,B07: real); stdcall; external archivo_dll;
   procedure InsercionBarras
       (var E5,     E4,D6,F5, D3,C8,H4,  H2,B5,E8,  F2,B7,G6, F4,D5,E6,
         E2,B8,H5,  H3,C4,D8,  G2,B6,F7,  D4,D7,G4: real); stdcall;  external archivo_dll;
   procedure PideBarras
       (var E5,     E4,D6,F5, D3,C8,H4,  H2,B5,E8,  F2,B7,G6, F4,D5,E6,
        E2,B8,H5,  H3,C4,D8,  G2,B6,F7,  D4,D7,G4: real); stdcall;  external archivo_dll;

   procedure UnPaso(var T,DT: real); stdcall; external archivo_dll;
   procedure Imprimir_Salida(var TTime,DT: real); stdcall; external archivo_dll;
   procedure Imprimir_Tiempos_CPU(var T1,T2,T3,T4: real); stdcall; external archivo_dll;
   procedure pedir_parametros_puntuales
     (var PotTotal0, PotResidual0, RO0,  KE0, DerivadaLogaritmica0, FForma0, PMax0: real;
      var CanalPMax0, TrozoPMax0: integer;
      var PotMax0: real; var CanalPotMax0: integer;
      var ConcXEMedia0, ConcXenon00: real;
      var ConcSMMedia0, ConcSamario00: real;
      var NumIter: Integer; var EPS,EPSKE: real) ; stdcall; external archivo_dll;

   procedure pasar_parametros_TH(var TC, TR, DR: CoreMatrix) ; stdcall; external archivo_dll;
   procedure pasar_BORO(var Boro: CoreMatrix) ; stdcall; external archivo_dll;
   procedure pedir_distribucion_potencias(var P: CoreMatrix; var PC: ChannelType) ; stdcall; external archivo_dll;
   procedure pedir_flujos(var F:  CoreFluxMatrix) ; stdcall; external archivo_dll;
   procedure Transferir(var TTT: real; var grabar: integer; var buffer: MemType); stdcall; external archivo_dll;
   procedure inicializar_DLL; Stdcall; external archivo_dll;
   procedure Cerrar_Salida; stdcall;   external archivo_dll;
{$ENDIF}

   (*******************************************************************************)
   (*******************************************************************************)
   (*******************************************************************************)
   (*******************************************************************************)
   var
     TT0Init: real;
     TFinal: real;
     TTime: real;
     DeltaT: real;
     TCPU_no_estacionario: real;
     TCPU_por_ciclo: real;
     TCPU_max_por_ciclo: real;
     TCPU_total: real;
     grabar, leer: boolean;
     grabando: integer;
     buffer: MemType;
     SF: file;



   (*******************************************************************************)
    function hacer_strr(a: single; n,d: integer) : string;
    var
      s: string;

    begin
      str (a:n:d,s);
      hacer_strr := s;
      end;   (*  hacer_strr   *)

   (*******************************************************************************)
    procedure PosicionarBarras (t: real);
    var
      B: real;
      B01,B02,B09,B11,B13,B03,B08,B10,B12,B07: real;

    begin
      if t <= 2.0 then
        B := 0.0
      else if (t <= 10.0) and (t>2.0) then
        B := 20.0*(t-2)/(10-2)
      else
        B := 20.0;

      B01 := 70.0;
      B02 := 110.0;
      B09 := B;
      B11 := 0.0;
      B13 := 0.0;
      B03 := 0.0;
      B08 := 0.0;
      B10 := 0.0;
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

    begin
      TTInit := Time;
      NIT := 0;
      DuracionMax := 0.0;

      while TTime < 0.999999*TFinal do begin
        TTime := TTime + DeltaT;

        TTInitCiclo := Time;
   //     PosicionarBarras (0.0);
        PosicionarBarras (TTime);

        UnPaso (TTime,DeltaT);

        DuracionCiclo := Time - TTInitCiclo;
        DuracionMax := max(DuracionCiclo,DuracionMax);

        imprimir_salida(TTime,DeltaT);

        inc(NIT);
        end;

      imprimir_salida(TTime,DeltaT);
      NIT := max(NIT,1);

      TCPU_no_estacionario := (Time-TTInit)*86400.0;
      TCPU_por_ciclo := (Time-TTInit)*86400.0/NIT;
      TCPU_max_por_ciclo := DuracionMax*86400.0;

      end;

   (*******************************************************************************)
    procedure EstadoEstacionario;
    var
      n, canal, trozo: integer;
      B01,B02,B09,B11,B13,B03,B08,B10,B12,B07: real;
      cero: real;

   {$IFDEF UseDLL}
      PPMBoro: CoreMatrix;
   {$ENDIF}

    begin

      B01 := 70.0;    B02 :=110 ; B09 := 0.0; B11 := 0.0;   B13 := 0.0;
      B03 := 0.0;   B08 := 0;   B10 := 70.0 ;   B12 :=00.0;   B07 := 0.0;

      cero := 0.0;

      for canal:=1 to NCANALES do
      for trozo:=1 to NTROZOS do
        PPMBoro[canal,trozo] := 0.0;

   {$IFDEF UseDLL}
      Pasar_BORO(PPMBoro);
   {$ENDIF}

      InsercionBancos  (B01,B02,B09,B11,B13,B03,B08,B10,B12,B07);


      for n:=1 to 10 do begin
        UnPaso(cero,cero);

        end;
        imprimir_salida(cero,cero);

      end;

   (*******************************************************************************)
   (*******************************************************************************)
   (*******************************************************************************)
   (*******************************************************************************)
    begin
     TTime := 0.0;
     TT0Init := Time;
     DeltaT := 0.1;
     TCPU_no_estacionario := 0.0;
     TCPU_por_ciclo := 0.0;
     TCPU_max_por_ciclo := 0.0;
     TFinal := 10.0;

     leer := false;
     grabar := false;
     inicializar_DLL;

     if leer then begin
       grabando := 0;
       AssignFile(SF, 'saved.xxx');
       Reset(SF, 1);
       BlockRead(SF,buffer,SizeOf(Buffer));
       Close(SF);
       Transferir(TTime,  grabando, buffer);
       end

     else
       EstadoEstacionario;

     if TFINAL > 0.0 then
       CalculoCinetico;

     if grabar then begin
       grabando := 1;
       Transferir(TTime,  grabando, buffer);
       AssignFile(SF, 'saved.xxx');
       Rewrite(SF, 1);	{ Record size = 1 }
       BlockWrite(SF,buffer,SizeOf(Buffer));
       Close(SF);
       end;

     TCPU_total := (Time-TT0Init)*86400.0;

     Imprimir_Tiempos_CPU(TCPU_no_estacionario, TCPU_por_ciclo,
                      TCPU_max_por_ciclo,TCPU_total);

     writeln ('Tiempo Total: ' +  hacer_strr (TCPU_total,1,3)) ;
     readln(TCPU_total);
     Cerrar_Salida;
     end.

