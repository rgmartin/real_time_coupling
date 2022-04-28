//{$DEFINE USARHEXA}
unit AlReactor;
interface
uses
  math, sysutils,

{$IFDEF USARHEXA}
   CAREMHEXA;
{$ELSE}
   CAREMTRI;
{$ENDIF}

type
   TString = array [0..256] of char;
   PString = ^TString;
//   MemType = array[0..2500000] of byte;

 procedure InsercionBancos  (B01,B02,B09,B11,B13,B03,B08,B10,B12,B07: real) ;
 procedure InsercionBarras
// BANCO 1   BANCO 2   BANCO 9   BANCO 11   BANCO 13   BANCO 3
     (E5,     E4,D6,F5, D3,C8,H4,  H2,B5,E8,  F2,B7,G6, F4,D5,E6,
// BANCO 8   BANCO 10   BANCO 12    BANCO 7
      E2,B8,H5,  H3,C4,D8,  G2,B6,F7,  D4,D7,G4: real);

 procedure PideBarras
// BANCO 1   BANCO 2   BANCO 9   BANCO 11   BANCO 13   BANCO 3
     (var E5,     E4,D6,F5, D3,C8,H4,  H2,B5,E8,  F2,B7,G6, F4,D5,E6,
// BANCO 8   BANCO 10   BANCO 12    BANCO 7
      E2,B8,H5,  H3,C4,D8,  G2,B6,F7,  D4,D7,G4: real);

 procedure UnPaso(T,DT: real);
 procedure Imprimir_Salida(TTime,DT: real);
 procedure Imprimir_Tiempos_CPU(T1,T2,T3,T4: real);
 procedure pasar_parametros_TH(var TC, TR, DR: CoreMatrix);
 procedure pedir_parametros_TH(var TC, TR, DR: CoreMatrix);
 procedure Transferir(var TTT: real; var grabar: integer; var buffer: MemType);
 procedure pedir_distribucion_potencias(var P: CoreMatrix; var PC: ChannelType) ;
 procedure pasar_BORO(var Boro: CoreMatrix) ;
 procedure pedir_parametros_puntuales
  (var PotTotal0, PotResidual0, RO0, KE0, DerivadaLogaritmica0, FForma0, PMax0: real;
    var CanalPMax0, TrozoPMax0: integer;
    var PotMax0: real; var CanalPotMax0: integer;
    var ConcXEMedia0, ConcXenon00: real;
    var ConcSMMedia0, ConcSamario00: real;
    var NumIter: Integer; var EPS,EPSKE: real) ;

implementation

// ****************************************************************************
 procedure InsercionBancos  (B01,B02,B09,B11,B13,B03,B08,B10,B12,B07: real) ;
 begin
   BarrasPorBanco  (B01,B02,B09,B11,B13,B03,B08,B10,B12,B07) ;
   end;

// ****************************************************************************
 procedure InsercionBarras
     (E5,     E4,D6,F5, D3,C8,H4,  H2,B5,E8,  F2,B7,G6, F4,D5,E6,
      E2,B8,H5,  H3,C4,D8,  G2,B6,F7,  D4,D7,G4: real);
 begin
   IntroducirTodasLasBarras
     (E5,     E4,D6,F5, D3,C8,H4,  H2,B5,E8,  F2,B7,G6, F4,D5,E6,
      E2,B8,H5,  H3,C4,D8,  G2,B6,F7,  D4,D7,G4);

   end;

// ****************************************************************************
 procedure PideBarras
     (var E5,     E4,D6,F5, D3,C8,H4,  H2,B5,E8,  F2,B7,G6, F4,D5,E6,
      E2,B8,H5,  H3,C4,D8,  G2,B6,F7,  D4,D7,G4: real);

 begin
   PedirTodasLasBarras
     (E5,     E4,D6,F5, D3,C8,H4,  H2,B5,E8,  F2,B7,G6, F4,D5,E6,
      E2,B8,H5,  H3,C4,D8,  G2,B6,F7,  D4,D7,G4);

   end;

// ****************************************************************************
 procedure UnPaso(T,DT: real);
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
 procedure Imprimir_Salida(TTime,DT: real);
 begin
   salidas(TTime,DT);
   end;

// ****************************************************************************
 procedure Imprimir_Tiempos_CPU(T1,T2,T3,T4: real);
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

procedure pedir_parametros_puntuales
  (var PotTotal0, PotResidual0, RO0, KE0, DerivadaLogaritmica0, FForma0, PMax0: real;
    var CanalPMax0, TrozoPMax0: integer;
    var PotMax0: real; var CanalPotMax0: integer;
    var ConcXEMedia0, ConcXenon00: real;
    var ConcSMMedia0, ConcSamario00: real;
    var NumIter: Integer; var EPS,EPSKE: real) ;

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
(* procedure pedir_parametros_puntuales
  (var PotTotal0, PotResidual0, RO0, KE0, DerivadaLogaritmica0, FForma0, PMax0: real;
    var CanalPMax0, TrozoPMax0: integer;
    var PotMax0: real; var CanalPotMax0: integer;
    var ConcXEMedia0, ConcXenon00: real) ;

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
    end;
  *)

// ****************************************************************************
 procedure pasar_parametros_TH(var TC, TR, DR: CoreMatrix) ;
 begin
   TComb := TC;
   DensRefr := DR;
   TempRefr := TR;
   end;

// ****************************************************************************
 procedure pedir_parametros_TH(var TC, TR, DR: CoreMatrix) ;
 begin
   TC := TComb;
   DR := DensRefr;
   TR := TempRefr;
   end;

// ****************************************************************************
 procedure pasar_BORO(var Boro: CoreMatrix) ;
 begin
   PPMBoro := Boro;
   end;

// ****************************************************************************
 procedure pedir_distribucion_potencias(var P: CoreMatrix; var PC: ChannelType) ;
 begin
   p := POT;      // Distribución de potencia por canal y trozo
   PC := potcan;  // Potencias totales por canal
   end;

// ****************************************************************************
 procedure pedir_flujos(var F:  CoreFluxMatrix) ;
 begin
   F := flujo;    // Distribución de flujo para todos los grupos por canal y trozo
   end;

// ****************************************************************************
 procedure Transferir(var TTT: real; var grabar: integer; var buffer: MemType);
 begin
   TransferirEstado(TTT, (grabar = 1), buffer);
   end;

begin
{$IFDEF USARHEXA}
  ArchivoDeEntrada := 'ENTRADA_HEXA.TXT';
{$ELSE}
  ArchivoDeEntrada := 'ENTRADA.TXT';
{$ENDIF}
  ArchivoDeSalida := 'SALIDA.TXT';

  inicializar;

end.
