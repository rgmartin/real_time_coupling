library reactor_hexa;


uses
  sysutils,CAREMHEXA;

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
   salidas (TTime,DT);
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
// ****************************************************************************
procedure Cerrar_Salida; stdcall;
begin
  if not ArchivoSalidaAbierto then exit;

  writeln (F);
  writeln (F);
  writeln (F,'<FAB>Cierro la salida</FAB>');
  close (F);
  end;

exports
  InsercionBancos name 'InsercionBancos' ,
  InsercionBarras name 'InsercionBarras' ,
  UnPaso name 'UnPaso' ,
  Imprimir_Salida name 'Imprimir_Salida' ,
  Imprimir_Tiempos_CPU name 'Imprimir_Tiempos_CPU' ,
  pedir_parametros_puntuales name 'pedir_parametros_puntuales' ,
  pasar_parametros_TH name 'pasar_parametros_TH' ,
  pasar_BORO name 'pasar_BORO' ,
  pedir_distribucion_potencias name 'pedir_distribucion_potencias' ,
  pedir_flujos name 'pedir_flujos' ,
  PideBarras name 'PideBarras' ,
  Transferir name 'Transferir' ,
  pedir_parametros_TH name 'pedir_parametros_TH',
  inicializar_DLL name 'inicializar_DLL',
  Estacionario name 'Estacionario',
  Cerrar_Salida name 'Cerrar_Salida' ;

begin
  ArchivoDeEntrada := 'ENTRADA_HEXA.TXT';
  ArchivoDeSalida := 'SALIDA.TXT';
  ArchivoSalidaAbierto := false;
  inicializar;
end.


