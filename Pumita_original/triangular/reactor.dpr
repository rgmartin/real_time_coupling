library reactor;

{ Important note about DLL memory management: ShareMem must be the
  first unit in your library's USES clause AND your project's (select
  Project-View Source) USES clause if your DLL exports any procedures or
  functions that pass strings as parameters or function results. This
  applies to all strings passed to and from your DLL--even those that
  are nested in records and classes. ShareMem is the interface unit to
  the BORLNDMM.DLL shared memory manager, which must be deployed along
  with your DLL. To avoid using BORLNDMM.DLL, pass string information
  using PChar or ShortString parameters. }

uses
  sysutils,CAREMTRI;

type
  TString = array [0..256] of char;


// ****************************************************************************
 procedure InsercionBancos
     (var B01,B02,B09,B11,B13,B03,B08,B10,B12,B07: real); stdcall;
 begin
   BarrasPorBanco  (B01,B02,B09,B11,B13,B03,B08,B10,B12,B07) ;
   end;

// ****************************************************************************
 procedure InsercionBarras
     (var E5,     E4,D6,F5, D3,C8,H4,  H2,B5,E8,  F2,B7,G6, F4,D5,E6,
          E2,B8,H5,  H3,C4,D8,  G2,B6,F7,  D4,D7,G4: real); stdcall;
 begin
   IntroducirTodasLasBarras
     (E5,     E4,D6,F5, D3,C8,H4,  H2,B5,E8,  F2,B7,G6, F4,D5,E6,
          E2,B8,H5,  H3,C4,D8,  G2,B6,F7,  D4,D7,G4);

   end;

// ****************************************************************************
 procedure PideBarras
     (var E5,     E4,D6,F5, D3,C8,H4,  H2,B5,E8,  F2,B7,G6, F4,D5,E6,
      E2,B8,H5,  H3,C4,D8,  G2,B6,F7,  D4,D7,G4: real); stdcall;

 begin
   PedirTodasLasBarras
     (E5,     E4,D6,F5, D3,C8,H4,  H2,B5,E8,  F2,B7,G6, F4,D5,E6,
      E2,B8,H5,  H3,C4,D8,  G2,B6,F7,  D4,D7,G4);

   end;

// ****************************************************************************
 procedure UnPaso(var T,DT: real); stdcall;
 begin
   TTime := T;
   DeltaT := DT;
   if DT = 0.0 then
     StaticCalculation(100.0)

   else
     if MetodoAdiabatico then
       AdiabaticMethod(DeltaT,DeltaT)
     else
       DirectMethod(DeltaT,DeltaT);

   end;

// ****************************************************************************
 procedure Estacionario(var POT: real); stdcall;
 begin
   TTime := 0.0;
   StaticCalculation(POT);

   end;

// ****************************************************************************
 procedure Imprimir_Salida(var TTime,DT: real); stdcall;
 begin
   salidas(TTime,DT);
   end;

// ****************************************************************************
 procedure Imprimir_Tiempos_CPU(var T1,T2,T3,T4: real); stdcall;
 begin
   if not ArchivoSalidaAbierto then exit;

   writeln (F);
   writeln (F);
   writeln (F,'Tiempo (CPU) no estacionario  = ', T1:10:3,' seg');

   if T2 > 0.0 then begin
     writeln (F,'Tiempo (CPU) por ciclo        = ', T2:10:3,' seg');
     writeln (F,'Tiempo (CPU) máximo por ciclo = ', T3:10:3,' seg');
     end;

   writeln (F);
   write (F,'           Tiempo (CPU) total = ',T4:10:3,' seg');
   writeln (F);
   close (F);
   end;

// ****************************************************************************
 procedure pedir_parametros_puntuales
  (var PotTotal0, PotResidual0, RO0, KE0, DerivadaLogaritmica0, FForma0, PMax0: real;
    var CanalPMax0, TrozoPMax0: integer;
    var PotMax0: real; var CanalPotMax0: integer;
    var ConcXEMedia0, ConcXenon00: real;
    var ConcSMMedia0, ConcSamario00: real;
    var NumIter: Integer; var EPS,EPSKE: real) ; stdcall;

 begin
    PotTotal0 := PotTotal;          // Potencia térmica total
    PotResidual0 := PotResidual;    // Potencia de decaimiento
    RO0 := RO;                      // Reactividad
    KE0 := keffInit;
    DerivadaLogaritmica0 :=
       DerivadaLogaritmica;         // Derivada logarítmica en los cálculos cinéticos
    FForma0 := FForma;              // Factor de forma de la potencia
    PMax0 := PMax;                  // Máxima potencia específica
    CanalPMax0 := CanalPMax;        // Canal donde ocurre PMax
    TrozoPMax0 := TrozoPMax;        // Trozo donde ocurre PMax
    PotMax0 := PotMax;              // Potencia máxima por canal
    CanalPotMax0 := CanalPotMax;    // Canal de máxima potencia
    ConcXEMedia0 := ConcXEMedia;    // Concentración de Xenón media actual
    ConcXenon00 := ConcXenon0;      // Concentración de Xenón media inicial
    ConcSmMedia0 := ConcSmMedia;    // Concentración de Samario media actual
    ConcSamario00 := ConcSamario0;  // Concentración de Samario media inicial
    NumIter := NumIteraciones;      // Número de iteraciones en el cálculo espacial
    EPS := UltimaConvergencia;      // Último error al interrumpir el cálculo del flujo
    EPSKE := PrecisionKeff;         // Último error en KEFF al interrumpir el cálculo del flujo

    end;


// ****************************************************************************
 procedure pasar_parametros_TH(var TC, TR, DR: CoreMatrix) ; stdcall;
 begin
   TComb := TC;
   DensRefr := DR;
   TempRefr := TR;
   end;

// ****************************************************************************
 procedure pedir_parametros_TH(var TC, TR, DR: CoreMatrix) ; stdcall;
 begin
   TC := TComb;
   DR := DensRefr;
   TR := TempRefr;
   end;

// ****************************************************************************
 procedure pasar_BORO(var Boro: CoreMatrix) ; stdcall;
 begin
   PPMBoro := Boro;
   end;

// ****************************************************************************
 procedure pedir_distribucion_potencias(var P: CoreMatrix; var PC: ChannelType) ; stdcall;
 begin
   p := POT;      // Distribución de potencia por canal y trozo
   PC := potcan;  // Potencias totales por canal
   end;

// ****************************************************************************
 procedure pedir_flujos(var F:  CoreFluxMatrix) ; stdcall;
 begin
   F := flujo;    // Distribución de flujo para todos los grupos por canal y trozo
   end;

// ****************************************************************************
 procedure Transferir(var TTT: real; var grabar: integer; var buffer: MemType); stdcall;
 begin
   TransferirEstado(TTT, (grabar = 1), buffer);
   end;

// ****************************************************************************
 procedure inicializar_DLL; stdcall;
 begin
   inicializar;
   end;

// ****************************************************************************
// ****************************************************************************
exports
  InsercionBancos index 1 name 'InsercionBancos' ,
  InsercionBarras index 2 name 'InsercionBarras' ,
  UnPaso index 3 name 'UnPaso' ,
  Imprimir_Salida index 4 name 'Imprimir_Salida' ,
  Imprimir_Tiempos_CPU index 5 name 'Imprimir_Tiempos_CPU' ,
  pedir_parametros_puntuales index 6 name 'pedir_parametros_puntuales' ,
  pasar_parametros_TH index 7 name 'pasar_parametros_TH' ,
  pasar_BORO index 8 name 'pasar_BORO' ,
  pedir_distribucion_potencias index 9 name 'pedir_distribucion_potencias' ,
  pedir_flujos index 10 name 'pedir_flujos' ,
  PideBarras index 11 name 'PideBarras' ,
  Transferir index 12 name 'Transferir' ,
  pedir_parametros_TH index 13 name 'pedir_parametros_TH',
  inicializar_DLL index 14 name 'inicializar_DLL',
  Estacionario index 15 name 'Estacionario' ;

begin
  ArchivoDeEntrada := 'ENTRADA.TXT';
  ArchivoDeSalida := 'SALIDA.TXT';

  ArchivoSalidaAbierto := false;
  inicializar;
  end.
