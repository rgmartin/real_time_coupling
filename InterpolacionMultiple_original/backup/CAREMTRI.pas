(*$R-*)

unit CAREMTRI;
interface
uses Dialogs,
  math,
  interpol,
  sysutils;

const
  NTROZOS = 14;
  NCANALES = 61;
  ng = 5;

  RENUM: array[1..61] of integer =
                 (44,
                45,43,
               46,24,42,
              47,25,23,41,
            48,26,11,22,40,
             27,12,10,21,
            49,13,3, 9, 39,
              28,4, 2, 20,
            50,14,1, 8, 38,
             29, 5, 7,37,
            51,15, 6,19,61,
             30,16,18,36,
            52,31,17,35,60,
             53,32,34,59,
              54,33,58,
               55,57,
                56) ;

  NomCan: array[1..61] of string[2] =
   ('E5', 'F4', 'E4', 'D5', 'D6', 'E6', 'F5', 'G4', 'G3', 'F3',
    'E3', 'D4', 'C5', 'C6', 'C7', 'D7', 'E7', 'F6', 'G5', 'H3',
    'H2', 'G2', 'F2', 'E2', 'D3', 'C4', 'B5', 'B6', 'B7', 'B8',
    'C8', 'D8', 'E8', 'F7', 'G6', 'H5', 'H4', 'I3', 'I2', 'I1',
    'H1', 'G1', 'F1', 'E1', 'D2', 'C3', 'B4', 'A5', 'A6', 'A7',
    'A8', 'A9', 'B9', 'C9', 'D9', 'E9', 'F8', 'G7', 'H6', 'I5', 'I4');

  RENUM_INV: array[1..61] of integer =
  (31, 27, 22, 26, 35, 40, 36, 32, 23, 18,
   13, 17, 21, 30, 39, 44, 49, 45, 41, 28,
   19, 14, 9, 5, 8, 12, 16, 25, 34, 43,
   48, 53, 57, 54, 50, 46, 37, 33, 24, 15,
   10, 6, 3, 1, 2, 4, 7, 11, 20, 29,
   38, 47, 52, 56, 59, 61, 60, 58, 55, 51,42);

type
  GroupType  = array[1..ng] of real;
  ScattType  = array[1..ng,1..ng] of real;
  CoreMatrix = array[1..NCANALES,1..NTROZOS] of real;
  CoreFluxMatrix = array[1..NCANALES,1..NTROZOS] of GroupType;
  CoreScattMatrix = array[1..NCANALES,1..NTROZOS] of ScattType;
  ChannelType = array[1..NCANALES] of real;
  IntCoreMatrix = array[1..NCANALES,1..NTROZOS] of integer;
  MemType = array[0..2500000] of byte;

var
  F: text;
  ArchivoSalidaAbierto: boolean;           // Si fue abierto el archivo de salida F
                                           // Si otro programa lo abre debe dar a
 	                                   // esta variable el valor TRUE
  ArchivoDeSalida: string;
  ArchivoDeEntrada: string;

  TTime,                                   // Instante en el que se está calculando un estado
  DeltaT,                                  // Paso de tiempo actual
  TFINAL: real;                            // Instante en el que se interrumpirá el cálculo
  PasoDeImpresion: real;                   // Paso para imprimir
  ProximaImpresion: real;                  // Instante para la próxima impresión

  NumIteracionesEstacionarias: integer;    // Número iteraciones para hallar el estado
                                           // estacionario
  NumIteraciones: integer;                 // Número total de iteraciones en el último
                                           // cálculo espacial del flujo
  UltimaConvergencia: real;                // Último error al interrumpir el cálculo del flujo
  PrecisionKeff:real;                      // Último error en KEFF al interrumpir el cálculo del flujo

  Cinetica: boolean;                       // Se simula un cálculo cinético espacial
  MetodoAdiabatico: boolean;
  Ciclo: boolean;                          // Se simula un ciclo de potencia

  GuardarEstado,                           // Guardar el último estado en el
                                           // archivo "SaveFileName"
  LeerEstado: boolean;                     // Caso de continuación en el que se lee
                                           // el primer estado de "RetrieveFileName"
  SaveFileName: string;                    // Archivo de para guardar el último estado
  RetrieveFileName: string;                // Archivo de para leer el último estado

  PotTotal,             // Potencia térmica total
  PotResidual,          // Potencia de decaimiento
  PotInit,              // Potencia inicial
  RO,                   // Reactividad
  KeffInit: real;       // K efectivo del primer estado estacionario.
  DerivadaLogaritmica,  // Derivada logarítmica en los cálculos cinéticos.
  FForma,               // Factor de forma de la potencia

  PMax: real;           // Máxima potencia específica

  CanalPMax, TrozoPMax: integer; // Canal y trozo donde ocurre la
                                 // máxima potencia específica
  PotMax: real;                  // Potencia máxima por canal
  CanalPotMax: integer;          // Canal de máxima potencia

  ConcXEMedia: real;   // Concentración de Xenón media actual
  ConcXenon0: real;    // Concentración de Xenón media inicial

  ConcSmMedia: real;   // Concentración de Samario media actual
  ConcSamario0: real;  // Concentración de Samario media inicial

  TComb: CoreMatrix;     // Distribución de la temperatura del combustible
  DensRefr: CoreMatrix;  // Distribución de la densidad del refrigerante
  TempRefr: CoreMatrix;  // Distribución de la temperatura del refrigerante
  PPMBoro: CoreMatrix;   // Distribución de la concentración de boro en ppm
  POT: CoreMatrix;       // Distribución de potencia por canal y trozo
  potcan: ChannelType;         // Potencias totales por canal

  flujo: CoreFluxMatrix;

 procedure inicializar;
 procedure StaticCalculation(POWER: real);
 procedure XenonStaticCalculation(POWER, DTXenon: real);
 procedure DirectMethod(DT, DTXenon: real);
 procedure AdiabaticMethod(DT,DTXenon: real);
 procedure salidas (TTime, DT: real);
 procedure TransferirEstado(var TTT: real; grabar: boolean; var mem: MemType);
 procedure BarrasPorBanco  (B01,B02,B09,B11,B13,B03,B08,B10,B12,B07: real) ;
 procedure IntroducirTodasLasBarras
// BANCO 1   BANCO 2   BANCO 9   BANCO 11   BANCO 13   BANCO 3
     (E5,     E4,D6,F5, D3,C8,H4,  H2,B5,E8,  F2,B7,G6, F4,D5,E6,
// BANCO 8   BANCO 10   BANCO 12    BANCO 7
      E2,B8,H5,  H3,C4,D8,  G2,B6,F7,  D4,D7,G4: real);
procedure PedirTodasLasBarras
// BANCO 1   BANCO 2   BANCO 9   BANCO 11   BANCO 13   BANCO 3
   (var E5,     E4,D6,F5, D3,C8,H4,  H2,B5,E8,  F2,B7,G6, F4,D5,E6,
// BANCO 8   BANCO 10   BANCO 12    BANCO 7
  E2,B8,H5,  H3,C4,D8,  G2,B6,F7,  D4,D7,G4: real);


(*********************************************************************************)
(*********************************************************************************)
(*********************************************************************************)
implementation

const
  NMaxTablas = 18;

Type
  TipoCalculoType = (estatico, adiabatico, directo);
  TipoSeccionesEficaces = array[1..40] of real;


var
 ArchivoIntMultiple: array[1..6] of string;
                                    // Archivo de interpolación múltiple.
                                    // Si no se da se supone interpolación en tablas.

 ConXenon,                          // Se utiliza reacoplamiento por xenón.
 ConSamario,                        // Se utiliza reacoplamiento por Samario.
 ConReacoplamientoTermohidraulico,  // Se supone reacoplamento por
	                            // temperatura de combustible,
	                            // densidad de refrigerante
	                            // temperatura del refrigerante
 ConBarras,                         // Se tiene en cuenta la inserción
	                            // de las barras de control
 ConQuemado,                        // Se tiene en cuenta la distribución de quemado
 ConBoro: boolean;                  // Se tiene en cuenta la distribución de boro

 variacion_lineal: boolean;         // Si se supone variación lineal de RO

 ImprimirIteraciones,      // Se muestra el proceso iterativo
 ImprimirNumIteraciones,   // Se da sólo el número de iteraciones
 ParaCompararConPUMA,      // Se da una lista de potencias por canal
	                   // para comparar con resultados de PUMA

 ImprimirDistribucion,     // Muestra la distribución de potencia por canal y trozo
 ImprimirPorCanal,         // Muestra la distribución de potencia por canal
 ImprimirTempDens:
       boolean;
 UNIDAD: real;

// Valores centrales de tabla
 ConcXE0: real;
 TComb0: real;
 DensRefr0: real;
 TempRefr0: real;

 PrecisionEstacionario, PrecisionCinetico: real;
 PuedeImprimir: boolean;
 FactPuma,KE: real;

// ********************************************** //
//   Variables y constantes de la interfase       //
// ********************************************** //

//   Inserciones de cada barra

//*** BANCO 1
  barra_E5: real;

//*** BANCO 2
  barra_E4, barra_D6, barra_F5: real ;

//*** BANCO 3
  barra_F4, barra_D5, barra_E6: real;

//*** BANCO 9
  barra_D3, barra_C8, barra_H4: real;

//*** BANCO 10
  barra_H3, barra_C4, barra_D8: real;

//*** BANCO 11
  barra_H2, barra_B5, barra_E8: real;

//*** BANCO 12
  barra_G2, barra_B6, barra_F7: real;

//*** BANCO 13
  barra_F2, barra_B7, barra_G6: real;

//*** BANCO 7
  barra_D4, barra_D7, barra_G4: real;

//*** BANCO 8
  barra_E2, barra_B8, barra_H5: real;

//*********************************************************************************
//*********************************************************************************

  KE0: real;            // K efectivo estático. En cinética espacial es
                        // el K efectivo del primer estado estacionario.
  BurnUp: CoreMatrix;    // Distribución de quemado
  ConcXE: CoreMatrix;    // Distribución de la concentración de xenón
  ConcI: CoreMatrix;     // Distribución de la concentración de iodo
  ConcPt: CoreMatrix;    // Distribución de la concentración de Prometeo
  ConcSm: CoreMatrix;    // Distribución de la concentración de Samario

  NumMaterial: IntCoreMatrix;       // Distribución de los números
                                    // de material por canal y trozo

  NumMaterialBarra: IntCoreMatrix;  // Número de material con barra
  FraccionBarra: CoreMatrix;        // Fracciones de insercion de barra

  PotRes: CoreMatrix;       // Distribución de la potencia de decaimiento
  PotInst: CoreMatrix;      // Distribución de la potencia instantánea (neutrónica)


const
  NX = 39;
  NY = 24;
  NRegionesReactor = 571;
  NTotalRegiones = 101;

//  NZ = 20;
//  NREFL1 = 4;
//  NREFL2 = 17;
//  NZ = 18;
//  NREFL1 = 3;
//  NREFL2 = 16;
  NZ = 16;
  NREFL1 = 2;
  NREFL2 = 15;
//  NZ = 14;
//  NREFL1 = 1;
//  NREFL2 = 14;

  NPREC = 6;
  NCOF = 7;
  LongTableLine0 = 326;
  NumMaxLineasTabla = 27;
  NPTOTAL = NZ*NX*NY*ng;


 type
   TINT11 = array[1..11] of integer;
   TLugar = array[1..NPTotal] of TINT11;
   PTLugar = ^TLugar;
   TREAL11 = array[1..11] of real;
   TValor = array[1..NPTotal] of TREAL11;
   PTValor = ^TValor;

 var
   lugar: PTLugar;
   valor: PTValor;

Type
  TableLine = array[1..LongTableLine0] of real;
  TableType = array[1..NumMaxLineasTabla] of TableLine;
  LatticeAxialLine = array[1..NZ] of real;
  AllLatticeMatrix = array[1..NZ,1..NY,1..NX] of GroupType;
  LatticeMatrix = array[1..NZ,1..NY,1..NX] of real;
  ScattMatrix = array[1..NZ,1..NY,1..NX] of ScattType;
  PrecursorType = array[1..NPREC] of real;
  BetaLamdaType = array[1..NPREC,1..ng] of real;
  PAllLatticeMatrix =  ^AllLatticeMatrix;
  OneDimension = array[1..NPTOTAL] of real;
  POneDimension = ^OneDimension;

const
  ALTURA = 8.0;
  DELTAX = 16.0/3.0;
  DELTAY = 16.0/3.0;
  DELTAZ = 10.0;
  APOTEMA = 16.0/6.0;
  LADO = 9.237604307034012;
  SUP_TRI = 36.95041722813605;
  SUP_LADO = 92.37604307034012;

const
   cof: array[1..7] of real =
   (3.733333E-03,6.000000E-04,3.666667E-05,4.166667E-06,1.500000E-07,
     1.333333E-08,8.666667E-10) ;

   cofexp: array[1..7] of real =
   (0.333333333,0.0333333333,3.333333E-03,3.333333E-04,3.333333E-05,
     3.333333E-06,3.333333E-07) ;

   espectro_de_fision: GroupType = (0.765153, 0.234484, 0.000363, 0.0, 0.0);

   BetaNuclear: PrecursorType = (4.17E-4,1.457E-3,1.339E-3,3.339E-3,8.97E-4,3.2E-4) ;
   Lambda: PrecursorType = (1.244E-2,3.063E-2,1.139E-1,3.079E-1,1.198,3.212) ;
   Fuente: real = 100.0;

 VEL: array[1..ng] of real = (1.910885E+09, 4.413117E+08, 5733925.69, 566564.116, 240619.248);

 CHIS: array[1..NPREC,1..ng] of real =
  ((0.131860723, 0.858153660, 9.986359E-03, 0, 0),
   (0.164914148, 0.825352477, 9.733892E-03, 0, 0),
   (0.102492868, 0.889305555, 8.202481E-03, 0, 0),
   (0.220731677, 0.773114707, 6.153698E-03, 0, 0),
   (0.181822667, 0.810229863, 7.947202E-03, 0, 0),
   (0.193439537, 0.793528244, 0.0130318594, 0, 0));


  reticula: array [1..147,1..4] of integer =

// *************************
// *****    NUCLEO     *****
// *************************

// *** FILAS 4,5

// ***  región 1
  ((19, 22,   4, 6),

// *** FILAS 5,6

// ***  region 2
   (16, 19,   5, 7),

// ***  region 3
   (22, 25,   5, 7),


// *** FILAS 6,7

// ***  region 4
   (13, 16,   6, 8),

// ***  region 5
   (19, 22,   6, 8),

// ***  region 6
   (25, 28,   6, 8),


// *** FILAS 7,8

// ***  region 7
   (10, 13,   7, 9),

// ***  region 8
   (16, 19,   7, 9),

// ***  region 9
   (22, 25,   7, 9),

// ***  region 10
   (28, 31,   7, 9),


// *** FILAS 8,9

// ***  region 11
   ( 7, 10,   8, 10),

// ***  region 12
   ( 13, 16,  8, 10),

// ***  region 13
   ( 19, 22,  8, 10),

// ***  region 14
   ( 25, 28,  8, 10),

// ***  region 15
   ( 31, 34,  8, 10),


//*** FILAS 9,10

// ***  region 16
   (10, 13,   9, 11),

// ***  region 17
   (16, 19,   9, 11),

// ***  region 18
   (22, 25,   9, 11),

// ***  region 19
   (28, 31,   9, 11),


//*** FILAS 10,11

// ***  region 20
   ( 7, 10,  10, 12),

// ***  region 21
   ( 13, 16, 10, 12),

// ***  region 22
   ( 19, 22, 10, 12),

// ***  region 23
   ( 25, 28, 10, 12),

// ***  region 24
   ( 31, 34, 10, 12),


// *** FILAS 11,12

// ***  region 25
   (10, 13,  11, 13),

// ***  region 26
   (16, 19,  11, 13),

// ***  region 27
   (22, 25,  11, 13),

// ***  region 28
   (28, 31,  11, 13),


// *** FILAS 12,13

// ***  region 29
   ( 7, 10,  12, 14),

// ***  region 30
   ( 13, 16, 12, 14),

// ***  region 31
   ( 19, 22, 12, 14),

// ***  region 32
   ( 25, 28, 12, 14),

// ***  region 33
   ( 31, 34, 12, 14),


// *** FILAS 13,14

// ***  region 34
   (10, 13,  13, 15),

// ***  region 35
   (16, 19,  13, 15),

// ***  region 36
   (22, 25,  13, 15),

// ***  region 37
   (28, 31,  13, 15),


// *** FILAS 14,15

// ***  region 38
   ( 7, 10,  14, 16),

// ***  region 39
   ( 13, 16, 14, 16),

// ***  region 40
   ( 19, 22, 14, 16),

// ***  region 41
   ( 25, 28, 14, 16),

// ***  region 42
   ( 31, 34, 14, 16),


// *** FILAS 15,16

// ***  region 43
   (10, 13,  15, 17),

// ***  region 44
   (16, 19,  15, 17),

// ***  region 45
   (22, 25,  15, 17),

// ***  region 46
   (28, 31,  15, 17),


// *** FILAS 16,17

// ***  region 47
   ( 7, 10,  16, 18),

// ***  region 48
   ( 13, 16, 16, 18),

// ***  region 49
   ( 19, 22, 16, 18),

// ***  region 50
   ( 25, 28, 16, 18),

// ***  region 51
   ( 31, 34, 16, 18),


// *** FILAS 17,18

// ***  region 52
   (10, 13,  17, 19),

// ***  region 53
   (16, 19,  17, 19),

// ***  region 54
   (22, 25,  17, 19),

// ***  region 55
   (28, 31,  17, 19),


// *** FILAS 18,19

// ***  region 56
   ( 13, 16, 18, 20),

// ***  region 57
   ( 19, 22, 18, 20),

// ***  region 58
   ( 25, 28, 18, 20),



// *** FILAS 19,20

// ***  region 59
   (16, 19,  19, 21),

// ***  region 60
   (22, 25,  19, 21),


// *** FILAS 20,21

// ***  region 61
   (19, 22,  20, 22),


// *************************
// *****   REFLECTOR   *****
// *************************

// *** REGION 62
( 17 , 24 , 2 , 3 ) ,

// *** REGION 63
( 12 , 29 , 3 , 4 ) ,

// *** REGION 64
( 9 , 19 , 4 , 5 ) ,

// *** REGION 65
( 22 , 32 , 4 , 5 ) ,

// *** REGION 66
( 8 , 16 , 5 , 6 ) ,

// *** REGION 67
( 25 , 33 , 5 , 6 ) ,

// *** REGION 68
( 5 , 13 , 6 , 7 ) ,

// *** REGION 69
( 28 , 36 , 6 , 7 ) ,

// *** REGION 70
( 4 , 10 , 7 , 8 ) ,

// *** REGION 71
( 31 , 37 , 7 , 8 ) ,

// *** REGION 72
( 3 , 7 , 8 , 9 ) ,

// *** REGION 73
( 34 , 38 , 8 , 9 ) ,

// *** REGION 74
( 3 , 7 , 9 , 10 ) ,

// *** REGION 75
( 34 , 38 , 9 , 10 ) ,

// *** REGION 76
( 3 , 7 , 10 , 11 ) ,

// *** REGION 77
( 34 , 38 , 10 , 11 ) ,

// *** REGION 78
( 2 , 7 , 11 , 12 ) ,

// *** REGION 79
( 34 , 39 , 11 , 12 ) ,

// *** REGION 80
( 2 , 7 , 12 , 13 ) ,

// *** REGION 81
( 34 , 39 , 12 , 13 ) ,

// *** REGION 82
( 2 , 7 , 13 , 14 ) ,

// *** REGION 83
( 34 , 39 , 13 , 14 ) ,

// *** REGION 84
( 2 , 7 , 14 , 15 ) ,

// *** REGION 85
( 34 , 39 , 14 , 15 ) ,

// *** REGION 86
( 3 , 7 , 15 , 16 ) ,

// *** REGION 87
( 34 , 38 , 15 , 16 ) ,

// *** REGION 88
( 3 , 7 , 16 , 17 ) ,

// *** REGION 89
( 34 , 38 , 16 , 17 ) ,

// *** REGION 90
( 3 , 7 , 17 , 18 ) ,

// *** REGION 91
( 34 , 38 , 17 , 18 ) ,

// *** REGION 92
( 4 , 10 , 18 , 19 ) ,

// *** REGION 93
( 31 , 37 , 18 , 19 ) ,

// *** REGION 94
( 5 , 13 , 19 , 20 ) ,

// *** REGION 95
( 28 , 36 , 19 , 20 ) ,

// *** REGION 96
( 8 , 16 , 20 , 21 ) ,

// *** REGION 97
( 25 , 33 , 20 , 21 ) ,

// *** REGION 98
( 9 , 19 , 21 , 22 ) ,

// *** REGION 99
( 22 , 32 , 21 , 22 ) ,

// *** REGION 100
( 12 , 29 , 22 , 23 ) ,

// *** REGION 101
( 17 , 24 , 23 , 24 ) ,

// CONTORNO EXTERNO

// *** REGION 102
( 1 , 40 , 1 , 2 ) ,

// *** REGION 103
( 1 , 17 , 2 , 3 ) ,

// *** REGION 104
( 24 , 40 , 2 , 3 ) ,

// *** REGION 105
( 1 , 12 , 3 , 4 ) ,

// *** REGION 106
( 29 , 40 , 3 , 4 ) ,

// *** REGION 107
( 1 , 9 , 4 , 5 ) ,

// *** REGION 108
( 32 , 40 , 4 , 5 ) ,

// *** REGION 109
( 1 , 8 , 5 , 6 ) ,

// *** REGION 110
( 33 , 40 , 5 , 6 ) ,

// *** REGION 111
( 1 , 5 , 6 , 7 ) ,

// *** REGION 112
( 36 , 40 , 6 , 7 ) ,

// *** REGION 113
( 1 , 4 , 7 , 8 ) ,

// *** REGION 114
( 37 , 40 , 7 , 8 ) ,

// *** REGION 115
( 1 , 3 , 8 , 9 ) ,

// *** REGION 116
( 38 , 40 , 8 , 9 ) ,

// *** REGION 117
( 1 , 3 , 9 , 10 ) ,

// *** REGION 118
( 38 , 40 , 9 , 10 ) ,

// *** REGION 119
( 1 , 3 , 10 , 11 ) ,

// *** REGION 120
( 38 , 40 , 10 , 11 ) ,

// *** REGION 121
( 1 , 2 , 11 , 12 ) ,

// *** REGION 122
( 39 , 40 , 11 , 12 ) ,

// *** REGION 123
( 1 , 2 , 12 , 13 ) ,

// *** REGION 124
( 39 , 40 , 12 , 13 ) ,

// *** REGION 125
( 1 , 2 , 13 , 14 ) ,

// *** REGION 126
( 39 , 40 , 13 , 14 ) ,

// *** REGION 127
( 1 , 2 , 14 , 15 ) ,

// *** REGION 128
( 39 , 40 , 14 , 15 ) ,

// *** REGION 129
( 1 , 3 , 15 , 16 ) ,

// *** REGION 130
( 38 , 40 , 15 , 16 ) ,

// *** REGION 131
( 1 , 3 , 16 , 17 ) ,

// *** REGION 132
( 38 , 40 , 16 , 17 ) ,

// *** REGION 133
( 1 , 3 , 17 , 18 ) ,

// *** REGION 134
( 38 , 40 , 17 , 18 ) ,

// *** REGION 135
( 1 , 4 , 18 , 19 ) ,

// *** REGION 136
( 37 , 40 , 18 , 19 ) ,

// *** REGION 137
( 1 , 5 , 19 , 20 ) ,

// *** REGION 138
( 36 , 40 , 19 , 20 ) ,

// *** REGION 139
( 1 , 8 , 20 , 21 ) ,

// *** REGION 140
( 33 , 40 , 20 , 21 ) ,

// *** REGION 141
( 1 , 9 , 21 , 22 ) ,

// *** REGION 142
( 32 , 40 , 21 , 22 ) ,

// *** REGION 143
( 1 , 12 , 22 , 23 ) ,

// *** REGION 144
( 29 , 40 , 22 , 23 ) ,

// *** REGION 145
( 1 , 17 , 23 , 24 ) ,

// *** REGION 146
( 24 , 40 , 23 , 24 ) ,

// *** REGION 147
(  1 , 40 , 24 , 25 ));


  XS_REFL_Inf0: TipoSeccionesEficaces =
  (2.8124E+00,  1.4410E+00,  7.5463E-01,  3.0428E-01,  1.7567E-01,  4.2124E-04,  8.0377E-02,  5.0538E-04,
   0.0000E+00,  0.0000E+00,  0.0000E+00,  7.9039E-06,  1.2073E-01,  1.1805E-05,  1.3915E-06,  0.0000E+00,
   0.0000E+00,  7.7335E-04,  1.0523E-01,  9.1564E-03,  0.0000E+00,  0.0000E+00,  5.1124E-04,  7.1569E-03,
   4.4781E-01,  0.0000E+00,  0.0000E+00,  6.1831E-06,  4.9096E-01,  1.4115E-02,  0.0000E+00,  0.0000E+00,
   0.0000E+00,  0.0000E+00,  0.0000E+00,  0.0000E+00,  0.0000E+00,  0.0000E+00,  0.0000E+00,  0.0000E+00);

  XS_REFL_Sup0: TipoSeccionesEficaces =
  (3.2555E+00,  1.6696E+00,  8.7387E-01,  3.5490E-01,  2.0361E-01,  3.6055E-04,  3.6326E-02,  8.2990E-05,
   2.9296E-10,  7.5586E-12,  0.0000E+00,  6.8191E-06,  9.9800E-03,  6.0324E-07,  6.6127E-08,  0.0000E+00,
   0.0000E+00,  6.6783E-04,  6.5427E-03,  3.8416E-04,  0.0000E+00,  0.0000E+00,  6.0190E-03,  6.1035E-03,
   4.6710E-02,  0.0000E+00,  0.0000E+00,  5.2449E-07,  1.2184E-01,  1.2071E-02,  0.0000E+00,  0.0000E+00,
   0.0000E+00,  0.0000E+00,  0.0000E+00,  0.0000E+00,  0.0000E+00,  0.0000E+00,  0.0000E+00,  0.0000E+00);

  XS_REFL_Rad0: TipoSeccionesEficaces =
  (1.7374E+00,  9.5902E-01,  3.6747E-01,  3.4981E-01,  2.9756E-01,  9.8896E-04,  3.6326E-02,  8.2990E-05,
   2.9296E-10,  7.5586E-12,  0.0000E+00,  7.9406E-04,  9.9800E-03,  6.0324E-07,  6.6127E-08,  0.0000E+00,
   0.0000E+00,  1.1324E-02,  6.5427E-03,  3.8416E-04,  0.0000E+00,  0.0000E+00,  6.0190E-03,  8.8651E-02,
   4.6710E-02,  0.0000E+00,  0.0000E+00,  5.2449E-07,  1.2184E-01,  1.9542E-01,  0.0000E+00,  0.0000E+00,
   0.0000E+00,  0.0000E+00,  0.0000E+00,  0.0000E+00,  0.0000E+00,  0.0000E+00,  0.0000E+00,  0.0000E+00);

(*******************************************************************************)
(*******************************************************************************)

type
  MatRespType = array[1..ng,1..ng] of real;

const
   alfaaxialInf0: MatRespType =
 ((-4.97211E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00),
  ( 1.84148E-01,-4.21518E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00),
  ( 4.33337E-02, 1.86478E-01,-2.68771E-01, 2.56083E-03, 9.60002E-04),
  ( 2.51561E-02, 5.14681E-02, 1.22382E-01,-2.55075E-01, 2.28995E-01),
  ( 1.08455E-02, 2.16496E-02, 4.17822E-02, 1.87520E-01,-2.70664E-01));

   alfaaxialSup0: MatRespType =
 ((-4.99999E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00),
  ( 1.82907E-01,-4.22387E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00),
  ( 4.14578E-02, 1.86035E-01,-2.67089E-01, 3.19515E-03, 1.21206E-03),
  ( 2.15478E-02, 4.89087E-02, 1.22441E-01,-2.43686E-01, 2.47177E-01),
  ( 7.95719E-03, 1.76834E-02, 3.66747E-02, 1.74610E-01,-2.89138E-01));

   Negro: MatRespType =
 ((-0.469220E00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00),
  ( 0.00000E+00,-0.469220E00, 0.00000E+00, 0.00000E+00, 0.00000E+00),
  ( 0.00000E+00, 0.00000E+00,-0.469220E00, 0.00000E+00, 0.00000E+00),
  ( 0.00000E+00, 0.00000E+00, 0.00000E+00,-0.469220E00, 0.00000E+00),
  ( 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-0.469220E00));



var
  COEF_axial_inf,  COEF_axial_sup: MatRespType;

  PotInstAnt: CoreMatrix;
  PAux: array [1..NCOF] of CoreMatrix;

  ConcPrecursor: array[1..6] of LatticeMatrix;
  FI, FIA, FIANT : AllLatticeMatrix;

  XS_REFL_Inf: TipoSeccionesEficaces;
  XS_REFL_Sup: TipoSeccionesEficaces;
  XS_REFL_Rad: TipoSeccionesEficaces;

  NZMIN, NZMAX: integer;

  conMatrizDeRespuesta: boolean;
  InterpolacionMultipleInicializada: boolean;
  NumFilas: array[1..18] of integer;

  PPorFision: real;
  MaxIteraciones: integer;

  BetaLambda: PrecursorType;
  espectro_dinamico: GroupType;
  espectro_estatico: GroupType;

  VOLTROZO: real;

  table: TableType;
  AllTables: array[1..NMaxTablas] of TableType;
  interpolacion_lineal: boolean;

  DIF, SIGMABS, NUSIGMAFIS,
  SIGMAFIS, NUSIGMAFIS0:  CoreFluxMatrix;
  SIGMABSXenon:  CoreFluxMatrix;
  SIGMABSSamario:  CoreFluxMatrix;
  SCATTERING: CoreScattMatrix;

  D, SIGTOT, NUSIGF, NUSIGF0,
  SIGF, LEFTSIDE, diag, diag0: AllLatticeMatrix;

  SIGSCATT: ScattMatrix;

  VolCeldilla: real;

  Coeff: array[1..ng,1..NZ,1..NY,1..NX,1..6] of real;

  FracPotResidual: real;

  BETA0: real;
  KE_ANT: real;
  LS0: array[1..ng] of real;
  BetaNuclearTotal: real;
  MeanLifeTime: real;

(***********************************************************************)
(***           para resolver el efecto cúspide                        **)
(***********************************************************************)

function cusp (I: real) : real;
  const
    a = 0.0812182741;
    b = 0.918781726;
begin
  I := 1.0 - I;
  Result := 0.5*a/b*(-1.0 + sqrt (1.0 + 4*b*I/(a*a)));
  Result := 1.0 - Result;
  end;


(***********************************************************************)
(***           invertir la matriz a por gauss-jordan                  **)
(***********************************************************************)

 procedure invert (var a: MatRespType);
 const
   n = 5;

 var
   i,j,k,ii,g,m,m1,mr,jaux: integer;
   e,t,t1: real;
   l: array [1..5] of integer;
   nlinea: integer;
   b: GroupType;

 begin

   for i:=1 to n do
     l[i] := i;

   e := 0.0;
   for i:=1 to n do begin
     t := 0.0;

     for j:=1 to n do
       t := t + abs(a[i,j]);

     if t > e then e := t;
     end;

   e:=1.0e-13*e;

   if e = 0.0 then begin
     writeln ('Matriz singular 1');
     halt(0);
     end;

   m := 1;
   for  k:=1 to n do begin
     t := 0.0;

     for i:=1 to n do begin
       t1 := abs(a[i,k]);

       if t1 >= t then begin
         t:=t1; m:=i;
         end;
       end;

     if t <= e then begin
       writeln ('Matriz singular 2');
       halt(0);
       end;

     for jaux:=1 to n do begin
       b[jaux] := a[k,jaux]; a[k,jaux] := a[m,jaux]; a[m,jaux] := b[jaux];
       end;

     m1:=l[m]; l[m]:=l[k]; l[k]:=m1;

     t:=1.0/a[k,k];
     for j:=1 to n do begin
       if j <> k then a[k,j] := a[k,j]*t;
       end;

     for i:=1 to n do begin
       if i <> k then
         for j:=1 to n do begin
           if j <>k then a[i,j] := a[i,j] - a[k,j]*a[i,k];
         end;
       end;

     for i:=1 to n do begin
       if i <> k then a[i,k] := -a[i,k]*t;
       end;

     a[k,k] := t;

     end;

   for k:=1 to n-1 do begin
     j := k;
     while not ((l[j] = k) or (j > n)) do
       j := j + 1;

     if j <> k then begin
       for g:=1 to n do begin
         b[g] := a[g,j]; a[g,j]:=a[g,k]; a[g,k]:= b[g];
         end;
       mr:=l[j]; l[j]:=l[k]; l[k]:=mr;
       end;

   end;

 end; // invert;

(*******************************************************************************)
 procedure hallar_contorno (D: GroupType; alfa: MatRespType; var beta: MatRespType);
 var
   ds,s: real;
   a, b, c: MatRespType;
   altura: real;
   i,j,k: integer;

 begin
   altura := DELTAZ;

   ds := 0.5*altura;
   s := SUP_TRI;

   for i:=1 to ng do
   for j:=1 to ng do
     if i = j then
       a[i,j] := 1.0 - ds/D[i]*alfa[i,i]
     else
       a[i,j] := -ds/D[i]*alfa[i,j];

   b := a;

   invert (b);

   for i:=1 to ng do
   for j:=1 to ng do
     c[i,j] := 0.0;

   for i:=1 to ng do
   for j:=1 to ng do
   for k:=1 to ng do
     c[i,j] := c[i,j] + a[i,k]*b[k,j];

   for i:=1 to ng do
     b[i,i] := b[i,i] - 1.0;

   for i:=1 to ng do
   for j:=1 to ng do
     b[i,j] := b[i,j]*s/ds*D[i];

   beta := b;
   end;

(*******************************************************************************)
 function hacer_strr(a: real; n,d: integer) : string;
 var
   s: string;

 begin
   str (a:n:d,s);
   hacer_strr := s;
   end;   (*  hacer_strr   *)

(*******************************************************************************)
function hacer_str(i,n: integer) : string;
  var
   s: string;

  begin
    str (i:n,s);
    hacer_str := s;

    end;   (*  hacer_str   *)

(*******************************************************************************)
procedure IntroducirTodasLasBarras
// BANCO 1   BANCO 2   BANCO 9   BANCO 11   BANCO 13   BANCO 3
   (E5,     E4,D6,F5, D3,C8,H4,  H2,B5,E8,  F2,B7,G6, F4,D5,E6,
// BANCO 8   BANCO 10   BANCO 12    BANCO 7
  E2,B8,H5,  H3,C4,D8,  G2,B6,F7,  D4,D7,G4: real);

begin

//*** BANCO 1
  barra_E5 := E5;

//*** BANCO 2
  barra_E4 := E4;  barra_D6 := D6;  barra_F5 := F5;

//*** BANCO 3
  barra_F4 := F4;  barra_D5 := D5;  barra_E6 := E6;

//*** BANCO 9
  barra_D3 := D3;  barra_C8 := C8;  barra_H4 := H4;

//*** BANCO 10
  barra_H3 := H3;  barra_C4 := C4;  barra_D8 := D8;

//*** BANCO 11
  barra_H2 := H2;  barra_B5 := B5;  barra_E8 := E8;

//*** BANCO 12
  barra_G2 := G2;  barra_B6 := B6;  barra_F7 := F7;

//*** BANCO 13
  barra_F2 := F2;  barra_B7 := B7;  barra_G6 := G6;

//*** BANCO 7
  barra_D4 := D4;  barra_D7 := D7;  barra_G4 := G4;

//*** BANCO 8
  barra_E2 := E2;  barra_B8 := B8;  barra_H5 := H5;

 end;

(*******************************************************************************)
procedure PedirTodasLasBarras
// BANCO 1   BANCO 2   BANCO 9   BANCO 11   BANCO 13   BANCO 3
   (var E5,     E4,D6,F5, D3,C8,H4,  H2,B5,E8,  F2,B7,G6, F4,D5,E6,
// BANCO 8   BANCO 10   BANCO 12    BANCO 7
  E2,B8,H5,  H3,C4,D8,  G2,B6,F7,  D4,D7,G4: real);

begin

//*** BANCO 1
  E5 := barra_E5;

//*** BANCO 2
  E4 := barra_E4;  D6 := barra_D6;  F5 := barra_F5;

//*** BANCO 3
  F4 := barra_F4;  D5 := barra_D5;  E6 := barra_E6;

//*** BANCO 9
  D3 := barra_D3;  C8 := barra_C8;  H4 := barra_H4;

//*** BANCO 10
  H3 := barra_H3;  C4 := barra_C4;  D8 := barra_D8;

//*** BANCO 11
  H2 := barra_H2;  B5 := barra_B5;  E8 := barra_E8;

//*** BANCO 12
  G2 := barra_G2;  B6 := barra_B6;  F7 := barra_F7;

//*** BANCO 13
  F2 := barra_F2;  B7 := barra_B7;  G6 := barra_G6;

//*** BANCO 7
  D4 := barra_D4;  D7 := barra_D7;  G4 := barra_G4;

//*** BANCO 8
  E2 := barra_E2;  B8 := barra_B8;  H5 := barra_H5;

 end;

(*******************************************************************************)
 procedure CompararConPuma (NumSegmento: integer);
 var
   iz,i,j,grupo, grupo1: integer;
   NomSalPuma, s: string;
   A,B,V: array[1..NY,1..NX] of real;

(* ------------------------------------------------------------------- *)
 procedure mostrar(tit: string);
 var
   i,j: integer;
   TodosNulos: boolean;

 begin
   writeln(F);
   writeln(F);
   writeln(F,tit);

   TodosNulos := true;

   for i:=1 to NY do
   for j:=1 to NX do
     if B[i,j] <> 0.0 then
       todosNulos := false;

   if TodosNulos then begin
     writeln(F);
     writeln(F,'Para ',tit,' Todos nulos');
     end

   else begin
     i := 1;
     while i <= NY do begin
       writeln(F);
       writeln(F,'Línea ',i);

       j := 1;
       while j <= NX do begin
         write(F,FloatToStrF(B[i,j],ffexponent,5,2),' ');
         if j mod 10 = 0 then writeln(F);
         inc(j);
         end;

       inc(i);
       end;

     end;

   end;

 (* ----------------------------------------------------- *)
 begin
   s := '';
   NomSalPuma := 'C:\ADAPUMA\SALIDA.';

   writeln (F);
   writeln (F);
   writeln (F,'Constantes de malla para el segmento ',NumSegmento);
   writeln (F);

   for i:=1 to NY do
   for j:=1 to NX do begin
     A[i,j] := 0.0;
     B[i,j] := 0.0;
     V[i,j] := 0.0
     end;

//   AssignFile(E, NomSalPuma);
//   Reset(E);

   for grupo:= 1 to ng do begin
     for i:=1 to NY do
     for j:=1 to NX do begin
       B[i,j] := coeff[grupo,NumSegmento,i,j,2];
       end;

     mostrar('AR0 ' + IntToStr(grupo));
     end;
// ----------------------------
   for grupo:= 1 to ng do begin

     for i:=1 to NY do
     for j:=1 to NX do begin
       B[i,j] := coeff[grupo,NumSegmento,i,j,3];
       end;
     mostrar('AT0' + IntToStr(grupo));
     end;

// ----------------------------
   for grupo:= 1 to ng do begin

     for i:=1 to NY do
     for j:=1 to NX do begin
       B[i,j] := coeff[grupo,NumSegmento,i,j,1];
       end;

     mostrar('AL0 ' + IntToStr(grupo));
     end;

// ----------------------------
   for grupo:= 1 to ng do begin

     for i:=1 to NY do
     for j:=1 to NX do begin
       B[i,j] := coeff[grupo,NumSegmento,i,j,4];
       end;

     mostrar('AB0 ' + IntToStr(grupo));
     end;

// ----------------------------
   for grupo:= 1 to ng do begin

     for i:=1 to NY do
     for j:=1 to NX do begin
       B[i,j] := coeff[grupo,NumSegmento,i,j,5];
       end;

     mostrar('AI0 ' + IntToStr(grupo));
     end;

// ----------------------------
   for grupo:= 1 to ng do begin

     for i:=1 to NY do
     for j:=1 to NX do begin
       B[i,j] := coeff[grupo,NumSegmento,i,j,6];
       end;

     mostrar('AS0 ' + IntToStr(grupo));
     end;


// ----------------------------
   for grupo:= 1 to ng do
   for grupo1:= 1 to ng do begin

     for i:=1 to NY do
     for j:=1 to NX do begin
       if grupo <> grupo1 then
         B[i,j] := SigSCATT[NumSegmento,i,j,grupo,grupo1]
       else
//         B[i,j] := Diag[NumSegmento,i,j,grupo];
         B[i,j] := SIGTOT[NumSegmento,i,j,grupo];
       end;

     mostrar('A5 ' + IntToStr(grupo) + ' ' + IntToStr(grupo1));
     end;

// ----------------------------
   for grupo:= 1 to ng do begin

     for i:=1 to NY do
     for j:=1 to NX do begin
       B[i,j] := NUSigF[NumSegmento,i,j,grupo];
       end;

     mostrar('A7 ' + IntToStr(grupo));
     end;


   Close(F);
   writeln ('Listo');
   Halt(0);

   end;

(*******************************************************************************)
  procedure MeshContants (TipoCalculo: TipoCalculoType; DT: real);
  const
    ALFA = 0.46922;
//    ALFA = 0.41300;

  var
    canal, trozo, grupo, grupo1, i, j, k, caso,
              N00, iz, i1, i2, j1, j2, prec: integer;
    LS: GroupType;
    ESCAPE_ABAJO, ESCAPE_ARRIBA: MatRespType;
    BetaLambdaDT: array[1..ng] of real;
    a1,a2,a3,x1,x2,x3,EX,VOL: real;
    New: boolean;
    rel1,rel2: real;
    DELTA: real;

  begin
    for grupo:=1 to ng do
     BetaLambdaDT[grupo] := 0.0;

    VOL := VOLCeldilla;

    if TipoCalculo = directo then begin
      for grupo:=1 to ng do
        LS0[grupo] := 1.0/(VEL[grupo]*DT);

      for grupo:=1 to ng do begin
        espectro_dinamico[grupo] :=
           (1.0-BetaNuclearTotal)*espectro_de_fision[grupo];
        for prec:=1 to nprec do
        espectro_dinamico[grupo] :=  espectro_dinamico[grupo]
           + lambda[prec]*chis[prec,grupo]*betanuclear[prec]*DT/(1.0 + lambda[prec]*DT);
        end;   
//      writeln(F);

//      for grupo:=1 to ng do
//        write(F,espectro_dinamico[grupo]:1:6,' ') ;

      end

    else
      for grupo:=1 to ng do
        LS0[grupo] := 0.0;

    for grupo:=1 to ng do
    for iz:=1 to NZ do
    for i:=1 to NY do
    for j:=1 to NX do
    if (D[iz,i,j,1] > 0.0)  and not ((iz < NZMIN) or (iz > NZMAX)) then begin
      SIGTOT[iz,i,j,grupo] :=  SIGTOT[iz,i,j,grupo]*VOL;
      NUSIGF[iz,i,j,grupo] :=  NUSIGF[iz,i,j,grupo]*VOL;
      for grupo1:=1 to ng do
        SIGSCATT[iz,i,j,grupo,grupo1] := SIGSCATT[iz,i,j,grupo,grupo1]*VOL;
      end;

    if ConMatrizDeRespuesta then
    for i:=1 to NY do
    for j:=1 to NX do
    if D[1,i,j,1] > 0.0 then begin
      hallar_contorno(D[1,i,j], alfaaxialInf0, ESCAPE_ABAJO);

      for grupo:=1 to ng do
      for grupo1:=1 to ng do
      if grupo = grupo1 then
        SIGTOT[1,i,j,grupo] := SIGTOT[1,i,j,grupo] - ESCAPE_ABAJO[grupo,grupo1]
      else
        SIGSCATT[1,i,j,grupo,grupo1] := SIGSCATT[1,i,j,grupo,grupo1]
                          +  ESCAPE_ABAJO[grupo,grupo1];

      end;

    if ConMatrizDeRespuesta then
    for i:=1 to NY do
    for j:=1 to NX do
    if D[NZ,i,j,1] > 0.0 then begin
      hallar_contorno(D[NZ,i,j], alfaaxialSup0, ESCAPE_ARRIBA);

      for grupo:=1 to ng do
      for grupo1:=1 to ng do
      if grupo = grupo1 then
        SIGTOT[NZ,i,j,grupo] := SIGTOT[NZ,i,j,grupo] - ESCAPE_ARRIBA[grupo,grupo1]
      else
        SIGSCATT[NZ,i,j,grupo,grupo1] := SIGSCATT[NZ,i,j,grupo,grupo1]
                          +  ESCAPE_ARRIBA[grupo,grupo1];
      end;

// CompararConPuma (5);

    for grupo:=1 to ng do
    for iz:=1 to NZ do
    for i:=1 to NY do
    for j:=1 to NX do

    if (D[iz,i,j,1] = 0.0)  or (iz < NZMIN) or (iz > NZMAX) then begin
      diag[iz,i,j,grupo] := 0.0;
      diag0[iz,i,j,grupo] := 0.0;

      for k:=1 to 6 do
        coeff[grupo,iz,i,j,k] := 0.0;

      end

    else begin

      diag[iz,i,j,grupo] := 0.0; //SIGTOT[grupo,iz,i,j];

      if (j > 1) and (D[iz,i,j-1,1] > 0.0) then begin
        EX := 2.0*SUP_LADO/(DELTAX/D[iz,i,j-1,grupo] + DELTAX/D[iz,i,j,grupo]);
        diag[iz,i,j,grupo] := diag[iz,i,j,grupo] + EX;
        coeff[grupo,iz,i,j,1] := EX;
        end

      else begin
        EX := SUP_LADO*ALFA/(1.0 + DELTAX*ALFA*0.5/D[iz,i,j,grupo]);
        coeff[grupo,iz,i,j,1] := EX;
        diag[iz,i,j,grupo] := diag[iz,i,j,grupo] + EX;
        end;


      if  (j < NX) and (D[iz,i,j+1,1] > 0.0) then begin
        EX := 2.0*SUP_LADO/(DELTAX/D[iz,i,j+1,grupo] + DELTAX/D[iz,i,j,grupo]);
        diag[iz,i,j,grupo] := diag[iz,i,j,grupo] + EX;
        coeff[grupo,iz,i,j,2] := EX;
        end

      else begin
        EX := (SUP_LADO*ALFA/(1.0 + DELTAX*0.5*ALFA/D[iz,i,j,grupo]));
        coeff[grupo,iz,i,j,2] := EX;
        diag[iz,i,j,grupo] := diag[iz,i,j,grupo] + EX;
        end;


      if (i + j) mod 2 = 0 then

      if (i > 1) and (D[iz,i-1,j,1] > 0.0) then begin
        EX := 2.0*SUP_LADO/(DELTAY/D[iz,i-1,j,grupo] + DELTAY/D[iz,i,j,grupo]);
        diag[iz,i,j,grupo] := diag[iz,i,j,grupo] + EX;
        coeff[grupo,iz,i,j,3] := EX;
        end

      else begin
        EX := SUP_LADO*ALFA/(1.0 + DELTAY*ALFA*0.5/D[iz,i,j,grupo]);
        coeff[grupo,iz,i,j,3] := EX;
        diag[iz,i,j,grupo] := diag[iz,i,j,grupo] + EX;
        end;

      if (i+ j) mod 2 = 1 then

      if (i < NY) and (D[iz,i+1,j,1] > 0.0) then begin
        EX := 2.0*SUP_LADO/(DELTAY/D[iz,i+1,j,grupo] + DELTAY/D[iz,i,j,grupo]);
        diag[iz,i,j,grupo] := diag[iz,i,j,grupo] + EX;
        coeff[grupo,iz,i,j,4] := EX;
        end

      else begin
        EX := SUP_LADO*ALFA/(1.0 + DELTAY*ALFA*0.5/D[iz,i,j,grupo]);
        coeff[grupo,iz,i,j,4] := EX;
        diag[iz,i,j,grupo] := diag[iz,i,j,grupo] + EX;
        end;

      if (iz > NZMIN) and (D[iz-1,i,j,1] > 0.0) then begin
        EX := 2.0*SUP_TRI/(DELTAZ/D[iz-1,i,j,grupo] + DELTAZ/D[iz,i,j,grupo]);
        diag[iz,i,j,grupo] := diag[iz,i,j,grupo] + EX;
        coeff[grupo,iz,i,j,5] := EX;
        end

      else if not ConMatrizDeRespuesta then begin
        EX := (SUP_TRI*ALFA/(1.0 + DELTAZ*ALFA*0.5/D[iz,i,j,grupo]));
        coeff[grupo,iz,i,j,5] := EX;
        diag[iz,i,j,grupo] := diag[iz,i,j,grupo] + EX;
        end;


      if (iz < NZMAX) and (D[iz+1,i,j,1] > 0.0) then begin
        EX := 2.0*SUP_TRI/(DELTAZ/D[iz+1,i,j,grupo] + DELTAZ/D[iz,i,j,grupo]);
        diag[iz,i,j,grupo] := diag[iz,i,j,grupo] + EX;
        coeff[grupo,iz,i,j,6] := EX;
        end

      else if not ConMatrizDeRespuesta then begin
        EX := (SUP_TRI*ALFA/(1.0 + DELTAZ*ALFA*0.5/D[iz,i,j,grupo]));
        coeff[grupo,iz,i,j,6] := EX;
        diag[iz,i,j,grupo] := diag[iz,i,j,grupo] + EX;
        end;

      end;

    diag0 := diag;

    for iz:=1 to NZ do
    for i:=1 to NY do
    for j:=1 to NX do
    for grupo:=1 to ng do
    if D[iz,i,j,grupo] > 0.0 then begin
      diag[iz,i,j,grupo] := diag0[iz,i,j,grupo] + SIGTOT[iz,i,j,grupo];
      end;

     for iz:=1 to NZ do
     for i:=1 to NY do
     for j:=1 to NX do
     if D[iz,i,j,1] > 0.0 then
     if TipoCalculo = directo then begin

       for grupo:=1 to ng do begin
         NUSIGF0[iz,i,j,grupo] :=NUSIGF0[iz,i,j,grupo]/KE0;
         NUSIGF[iz,i,j,grupo] :=NUSIGF[iz,i,j,grupo]/KE0;
         end;

       for grupo:=1 to ng do begin
         LS[grupo] := LS0[grupo]*FI[iz,i,j,grupo];

         for prec:=1 to NPREC do
           LS[grupo] := LS[grupo]
              + chis[prec,grupo]*Lambda[prec]/(1.0 + lambda[prec]*DT)*ConcPrecursor[prec,iz,i,j];

         LeftSide[iz,i,j,grupo] :=  (LS[grupo] + fuente)*VOL;

//         NUSIGF[iz,i,j,grupo] :=
//             NUSIGF[iz,i,j,grupo]*((1.0-BetaEfTotal) + BetaLambdaDT[grupo]);
         end


       end

     else if TipoCalculo = adiabatico then
       for grupo:=1 to ng do begin
         NUSIGF0[iz,i,j,grupo] :=NUSIGF0[iz,i,j,grupo]/KE0;
         NUSIGF[iz,i,j,grupo] :=NUSIGF[iz,i,j,grupo]/KE0;
         end;

    end;

(*******************************************************************************)
 procedure BarrasPorBanco  (B01,B02,B09,B11,B13,B03,B08,B10,B12,B07: real) ;
 begin
//*** BANCO 1
  barra_E5 := B01;

//*** BANCO 2
  barra_E4 := B02 ;
  barra_D6 := B02 ;
  barra_F5 := B02 ;

//*** BANCO 3
  barra_F4 := B03 ;
  barra_D5 := B03 ;
  barra_E6 := B03 ;

//*** BANCO 9
  barra_D3 := B09 ;
  barra_C8 := B09 ;
  barra_H4 := B09 ;

//*** BANCO 10
  barra_H3 := B10 ;
  barra_C4 := B10 ;
  barra_D8 := B10 ;

//*** BANCO 11
  barra_H2 := B11 ;
  barra_B5 := B11 ;
  barra_E8 := B11 ;

//*** BANCO 12
  barra_G2 := B12 ;
  barra_B6 := B12 ;
  barra_F7 := B12 ;

//*** BANCO 13
  barra_F2 := B13 ;
  barra_B7 := B13 ;
  barra_G6 := B13 ;

//*** BANCO 7
  barra_D4 := B07 ;
  barra_D7 := B07 ;
  barra_G4 := B07 ;

//*** BANCO 8
  barra_E2 := B08 ;
  barra_B8 := B08 ;
  barra_H5 := B08 ;


   end;

(*******************************************************************************)
  function InterpolateTable(BU: real; NMAT, n: integer;
                 var a1,a2,a3,x1,x2,x3: real; var New: boolean): real;
(* const
    XS: array[1..40] of real =
 (2.61698E+00, 1.36062E+00, 8.22178E-01, 4.40172E-01, 2.52215E-01,
  2.65897E-03, 6.43050E-02, 5.30847E-04, 1.24768E-11, 2.05254E-13,
  0.00000E+00, 1.43602E-03, 7.07928E-02, 4.15708E-06, 4.09713E-07,
  0.00000E+00, 0.00000E+00, 1.78107E-02, 5.22929E-02, 4.68823E-03,
  0.00000E+00, 0.00000E+00, 2.01360E-03, 5.36672E-02, 2.50017E-01,
  0.00000E+00, 0.00000E+00, 4.83469E-06, 3.51147E-01, 1.01826E-01,
  6.07205E-03, 5.86657E-04, 8.35188E-03, 6.70160E-02, 1.45868E-01,
  2.40071E-03, 2.52309E-04, 3.61949E-03, 2.90310E-02, 6.31894E-02);
  *)
  var
   i, LongTable: integer;

  begin
//    Result := XS[n-1];
//    new := false;
//    exit;


     LongTable := NumFilas[NMAT];
     i := 1;
     while (i <= LongTable) and (Table[i,1] < BU) do inc(i);

     if i > LongTable then Result := Table[LongTable,n]
     else if i = 1 then Result := Table[1,n]
     else begin
       if i = 1 then i := 2
       else if i = LongTable then i := LongTable-1;

       if interpolacion_lineal then begin

         if New then begin
           x1 := Table[i-1,1];
           x2 := Table[i,1];
           x3 := 0.0;
           a1 := (BU-x2)/(x1-x2);
           a2 := (BU-x1)/(x2-x1);
           a3 := 0.0;
           end;

         Result := Table[i-1,n]*a1 + Table[i,n]*a2;
         end

       else begin

         if New then begin
           x1 := Table[i-1,1];
           x2 := Table[i,1];
           x3 := Table[i+1,1];
           a1 := (BU-x2)*(BU-x3)/((x1-x2)*(x1-x3));
           a2 := (BU-x1)*(BU-x3)/((x2-x1)*(x2-x3));
           a3 := (BU-x1)*(BU-x2)/((x3-x1)*(x3-x2));
           end;

         Result := Table[i-1,n]*a1 + Table[i,n]*a2
                                  + Table[i+1,n]*a3;
         end;

       end;

     new := false;
     end;

(*******************************************************************************)
  procedure Imprimir_XS_Red (NX1,NX2,NY1,NY2,NZ1,NZ2: integer);
  var
    i,j,kz,grupo,grupo1: integer;
    ss: real;

  begin
    for kz:=NZ1 to NZ2 do begin
      writeln(F);
      writeln(F,'XS DE LA RED NIVEL ',kz);
      for i:=NY1 to NY2 do
      for j:=NX1 to NX2 do begin
        writeln(F);
        writeln(F,'fila ',i,'  col ',j);

        for grupo:=1 to ng do begin
           if D[kz,i,j,grupo] = 0.0 then
             ss := 0.0
           else
             ss := 1.0/(3.0*D[kz,i,j,grupo]);

           write(F,FloatToStrF(ss,ffexponent,5,2),' ');
           end;

        writeln(F);
        for grupo:=1 to ng do begin
          for grupo1:=1 to ng do
            if grupo1 = grupo then
              write(F,FloatToStrF(SIGTOT[kz,i,j,grupo],ffexponent,5,2),' ')
            else
              write(F,FloatToStrF(SigSCATT[kz,i,j,grupo,grupo1],ffexponent,5,2),' ');
          writeln(F);
          end;

        for grupo:=1 to ng do
          write(F,FloatToStrF(NUSIGF[kz,i,j,grupo],ffexponent,5,2),' ');

        writeln(F);
        for grupo:=1 to ng do
          write(F,FloatToStrF(SIGF[kz,i,j,grupo],ffexponent,5,2),' ');

        end;
      end;

   Close(F);
   writeln ('Listo');
   Halt(0);
   end;

(****************************************************************)
procedure materiales (canal: integer; tipo: string; VQ0: string;
              tipo_barra: integer; Ins: real);
var
  BARRA: boolean;
  trozo, num_dedos, NT, NumIns: integer;
  FracTrozo: real;
  VQ: string;
  Fraccionar, CASO: boolean;

begin
  NT := 1;
  NumIns := trunc(Ins) div 10 ;
  FracTrozo := frac(Ins/10.0);
  if FracTrozo > 0.0 then
    FracTrozo := cusp(FracTrozo);

  for trozo:=1 to 14 do begin
    BARRA := trozo >= 15 - NumIns ;

    if (trozo = 1) or (trozo >= 13) then VQ := ' 0  '
                                    else VQ := VQ0;
    Fraccionar := (FracTrozo > 0.0) and (trozo = 14 - NumIns);

    for CASO:=false to Fraccionar do begin
      if Fraccionar and caso then
        BARRA := true;

      if BARRA then
        num_dedos := tipo_barra
      else
        num_dedos := 0;

      if num_dedos = 0 then begin
        if tipo = '1.8%' then NT := 1;
        if tipo = '3.1%' then
          if VQ = ' 0  ' then NT := 3;
        if tipo = '3.1%' then
          if VQ = ' 6VQ' then NT := 4;
        end;

      if num_dedos = 12 then
        if tipo = '1.8%' then NT := 2;

      if num_dedos = 18 then begin
        if tipo = '3.1%' then
          if VQ = ' 0  ' then NT := 5;
        if tipo = '3.1%' then
          if VQ = ' 6VQ' then NT := 6 ;
        end;

//      if (trozo = 6) or (trozo = 11) then NT := NT + 6 ;

      if caso then begin
        NumMaterialBarra[canal,trozo] := NT;
        FraccionBarra[canal,trozo] := FracTrozo;
        end

      else begin
        NumMaterial[canal,trozo] := NT;
        NumMaterialBarra[canal,trozo] := 0;
        FraccionBarra[canal,trozo] := 0.0;
        end


      end;

    end;

  end;

(**************************************************************)
procedure asignar_materiales ;
const
  no_barra = 0;

begin
  materiales ( 1,'3.1%',' 6VQ', 0, no_barra) ;

  materiales ( 2,'3.1%',' 6VQ', 0, no_barra) ;
  materiales ( 3,'3.1%',' 6VQ', 0, no_barra) ;

  materiales ( 4,'3.1%',' 6VQ', 0, no_barra) ;
  materiales ( 5,'3.1%',' 6VQ',18, barra_E2) ;
  materiales ( 6,'3.1%',' 6VQ', 0, no_barra) ;

  materiales ( 7,'3.1%',' 6VQ', 0, no_barra) ;
  materiales ( 8,'3.1%',' 6VQ',18, barra_D3) ;
  materiales ( 9,'3.1%',' 6VQ',18, barra_F2) ;
  materiales (10,'3.1%',' 6VQ', 0, no_barra) ;

  materiales (11,'3.1%',' 6VQ', 0, no_barra) ;
  materiales (12,'3.1%',' 6VQ',18, barra_C4) ;
  materiales (13,'3.1%',' 6VQ', 0, no_barra) ;
  materiales (14,'3.1%',' 6VQ',18, barra_G2) ;
  materiales (15,'3.1%',' 6VQ', 0, no_barra) ;

  materiales (16,'3.1%',' 6VQ',18, barra_B5) ;
  materiales (17,'3.1%',' 6VQ',18, barra_D4) ;
  materiales (18,'3.1%',' 6VQ',0 , no_barra) ;
  materiales (19,'3.1%',' 6VQ',18, barra_H2) ;

  materiales (20,'3.1%',' 6VQ', 0, no_barra) ;
  materiales (21,'3.1%',' 6VQ', 0, no_barra) ;
  materiales (22,'3.1%',' 6VQ',18, barra_E4) ;
  materiales (23,'3.1%',' 6VQ', 0, no_barra) ;
  materiales (24,'3.1%',' 6VQ', 0, no_barra) ;

  materiales (25,'3.1%',' 6VQ',18, barra_B6) ;
  materiales (26,'3.1%',' 6VQ',18, barra_D5) ;
  materiales (27,'3.1%',' 6VQ',18, barra_F4) ;
  materiales (28,'3.1%',' 6VQ',18, barra_H3) ;

  materiales (29,'3.1%',' 6VQ', 0, no_barra) ;
  materiales (30,'3.1%',' 6VQ', 0, no_barra) ;
  materiales (31,'1.8%','   0',12, barra_E5) ;
  materiales (32,'3.1%',' 6VQ',18, barra_G4) ;
  materiales (33,'3.1%',' 6VQ', 0, no_barra) ;

  materiales (34,'3.1%',' 6VQ',18, barra_B7) ;
  materiales (35,'3.1%',' 6VQ',18, barra_D6) ;
  materiales (36,'3.1%',' 6VQ',18, barra_F5) ;
  materiales (37,'3.1%',' 6VQ',18, barra_H4) ;

  materiales (38,'3.1%',' 6VQ', 0, no_barra) ;
  materiales (39,'3.1%',' 6VQ', 0, no_barra) ;
  materiales (40,'3.1%',' 6VQ',18, barra_E6) ;
  materiales (41,'3.1%',' 6VQ', 0, no_barra) ;
  materiales (42,'3.1%',' 6VQ', 0, no_barra) ;

  materiales (43,'3.1%',' 6VQ',18, barra_B8) ;
  materiales (44,'3.1%',' 6VQ',18, barra_D7) ;
  materiales (45,'3.1%',' 6VQ', 0, no_barra) ;
  materiales (46,'3.1%',' 6VQ',18, barra_H5) ;

  materiales (47,'3.1%',' 6VQ', 0, no_barra) ;
  materiales (48,'3.1%',' 6VQ',18, barra_C8) ;
  materiales (49,'3.1%',' 6VQ', 0, no_barra) ;
  materiales (50,'3.1%',' 6VQ',18, barra_G6) ;
  materiales (51,'3.1%',' 6VQ', 0, no_barra) ;

  materiales (52,'3.1%',' 6VQ', 0, no_barra) ;
  materiales (53,'3.1%',' 6VQ',18, barra_D8) ;
  materiales (54,'3.1%',' 6VQ',18, barra_F7) ;
  materiales (55,'3.1%',' 6VQ', 0, no_barra) ;

  materiales (56,'3.1%',' 6VQ', 0, no_barra) ;
  materiales (57,'3.1%',' 6VQ',18, barra_E8) ;
  materiales (58,'3.1%',' 6VQ', 0, no_barra) ;

  materiales (59,'3.1%',' 6VQ', 0, no_barra) ;
  materiales (60,'3.1%',' 6VQ', 0, no_barra) ;

  materiales (61,'3.1%',' 6VQ', 0, no_barra) ;

  end;


(*******************************************************************************)
  procedure InterpolateXS ;
  const

   DBORO: array[1..40] of real =
  (-1.956000E-07, 3.594000E-07, 8.848000E-07, 1.182900E-06, 1.420900E-05,
    1.000000E-08,-2.600000E-07,-2.400000E-09, 4.000000E-17, 8.000000E-19,
    0           , 1.500000E-08, 7.400000E-07, 4.400000E-11, 5.100000E-12,
    0           , 0           , 8.200000E-07, 2.700000E-07, 2.500000E-08,
    0           , 0           , 2.130000E-07, 7.670000E-06,-6.100000E-06,
    0           , 0           , 1.000000E-12, 4.000000E-07, 1.741000E-05,
    2.800000E-08, 1.100000E-09, 3.100000E-08,-7.400000E-07,-3.000000E-07,
    9.000000E-09, 5.000000E-10, 1.300000E-08,-3.200000E-07,-1.100000E-07);



  var
    canal, trozo, i, j, k, N00, iz, i1, i2, j1, j2, pos, NMAT,
       grupo, grupo1: integer;
    a1,a2,a3,x1,x2,x3,Q, fact1, fact2, VAL: real;
    New, fraccionado, caso: boolean;
    DELTA, BORO, fraccion: real;
    P1, P2, P3, P4, P5, P6, P7: real;
    AllXS: TXS;
    ttt0: real;

(* ---------------------------------------------------------------------- *)
 procedure agregar_separador;
 const
   DSeparador: array[1..40] of real =
(*  (   2.558E-03,    6.061E-03,    7.018E-03,   -8.153E-03,   -2.392E-02,
     2.792000E-05,-6.611000E-04,-1.508400E-05,-2.606000E-09,-2.316000E-10,
     4.450000E-13, 3.830000E-05,-2.690300E-03,-1.625100E-07,-1.597300E-08,
         0       , 2.655300E-13, 2.412000E-04,-2.332200E-03,-2.163800E-04,
         0       ,     0       , 1.085300E-04, 2.014900E-03,-0.011765    ,
         0       ,     0       ,-1.858200E-07,   -0.011755 , 3.098600E-03,
    -4.117000E-05, 6.621000E-06, 8.980000E-05, 1.831300E-03, 1.768000E-03,
    -1.515000E-05, 2.616000E-06, 3.555000E-05, 7.028000E-04, 6.935000E-04);
  *)
   (1.279000E-03, 3.030500E-03, 3.509000E-03,-4.076500E-03,-0.01196     ,
    1.396000E-05,-3.305500E-04,-7.542000E-06,-1.303000E-09,-1.158000E-10,
    2.225000E-13, 1.915000E-05,-1.345150E-03,-8.125500E-08,-7.986500E-09,
    0           , 1.327650E-13, 1.206000E-04,-1.166100E-03,-1.081900E-04,
    0           , 0           , 5.426500E-05, 1.007450E-03,-5.882500E-03,
    0           , 0           ,-9.291000E-08,-5.877500E-03, 1.549300E-03,
   -2.058500E-05, 3.310500E-06, 4.490000E-05, 9.156500E-04, 8.840000E-04,
   -7.575000E-06, 1.308000E-06, 1.777500E-05, 3.514000E-04, 3.467500E-04);

 var
   grupo, grupo1: integer;

 begin
   pos := 1;

   for grupo:=1 to ng do begin
     DIF[canal,trozo,grupo] := DIF[canal,trozo,grupo] + DSeparador[pos];
     inc(pos);
     end;

   for grupo:=1 to ng do for grupo1:=1 to ng do begin
     if grupo = grupo1 then
       SIGMABS[canal,trozo,grupo] :=
         SIGMABS[canal,trozo,grupo] + DSeparador[pos]
    else
      SCATTERING[canal,trozo,grupo,grupo1] :=
         SCATTERING[canal,trozo,grupo,grupo1] + DSeparador[pos];
    inc(pos);
    end;

  for grupo:=1 to ng do begin
    NUSIGMAFIS[canal,trozo,grupo] := NUSIGMAFIS[canal,trozo,grupo] +
         DSeparador[pos];
    inc(pos);
    end;

  for grupo:=1 to ng do begin
    SIGMAFIS[canal,trozo,grupo] := SIGMAFIS[canal,trozo,grupo] +
         DSeparador[pos];
    inc(pos);
    end;

  end; (* agregar_separador *)

(* ---------------------------------------------------------------------- *)
 procedure agregar_boro;
 var
   grupo, grupo1: integer;

 begin
   DELTA := PPMBoro[canal,trozo];
   pos := 1;

   for grupo:=1 to ng do begin
     DIF[canal,trozo,grupo] := DIF[canal,trozo,grupo] + DELTA*DBORO[pos];
     inc(pos);
     end;

   for grupo:=1 to ng do for grupo1:=1 to ng do begin
     if grupo = grupo1 then
       SIGMABS[canal,trozo,grupo] :=
         SIGMABS[canal,trozo,grupo] + DELTA*DBORO[pos]
    else
      SCATTERING[canal,trozo,grupo,grupo1] :=
         SCATTERING[canal,trozo,grupo,grupo1] + DELTA*DBORO[pos];
    inc(pos);
    end;

  for grupo:=1 to ng do begin
    NUSIGMAFIS[canal,trozo,grupo] := NUSIGMAFIS[canal,trozo,grupo] +
         DELTA*DBORO[pos];
    inc(pos);
    end;

  for grupo:=1 to ng do begin
    SIGMAFIS[canal,trozo,grupo] := SIGMAFIS[canal,trozo,grupo] +
         DELTA*DBORO[pos];
    inc(pos);
    end;

  end; (* agregar_boro *)

(* ---------------------------------------------------------------------- *)
 procedure agregar (N00: integer; DELTA: real);
 var
   grupo, grupo1, pos: integer;
 begin
   pos := N00+1;
   for grupo:=1 to ng do begin
     DIF[canal,trozo,grupo] := DIF[canal,trozo,grupo] +
       DELTA*InterpolateTable(Q,NMAT,pos,a1,a2,a3,x1,x2,x3,New);
     inc(pos);
     end;

   for grupo:=1 to ng do for grupo1:=1 to ng do begin
     if grupo = grupo1 then
       SIGMABS[canal,trozo,grupo] := SIGMABS[canal,trozo,grupo] +
         DELTA*InterpolateTable(Q,NMAT,pos,a1,a2,a3,x1,x2,x3,New)
     else
       SCATTERING[canal,trozo,grupo,grupo1] := SCATTERING[canal,trozo,grupo,grupo1] +
           DELTA*InterpolateTable(Q,NMAT,pos,a1,a2,a3,x1,x2,x3,New);
     inc(pos);
     end;

   for grupo:=1 to ng do begin
     NUSIGMAFIS[canal,trozo,grupo] := NUSIGMAFIS[canal,trozo,grupo] +
              DELTA*InterpolateTable(Q,NMAT,pos,a1,a2,a3,x1,x2,x3,New);
     inc(pos);
     end;

   for grupo:=1 to ng do begin
     SIGMAFIS[canal,trozo,grupo] := SIGMAFIS[canal,trozo,grupo] +
              DELTA*InterpolateTable(Q,NMAT,pos,a1,a2,a3,x1,x2,x3,New);
     inc(pos);
     end;

   end;  // agregar

(* ---------------------------------------------------------------------- *)
 procedure hallar_SIGMABS_Xenon (N00: integer);
 var
   grupo, grupo1, pos: integer;
 begin
   pos := N00+ng+1;
   new := true;

   for grupo:=1 to ng do for grupo1:=1 to ng do begin
     if grupo = grupo1 then begin
       SIGMABSXenon[canal,trozo,grupo] :=
         InterpolateTable(Q,NMAT,pos,a1,a2,a3,x1,x2,x3,New);
       end;
     inc(pos);
     end;

   end;  // hallar_SIGMABS_Xenon

(* ---------------------------------------------------------------------- *)
 procedure agregar_Xenon (N00: integer; DELTA: real);
 var
   grupo, grupo1, pos: integer;
 begin
   pos := N00+1;
   for grupo:=1 to ng do begin
     DIF[canal,trozo,grupo] := DIF[canal,trozo,grupo] +
       DELTA*InterpolateTable(Q,NMAT,pos,a1,a2,a3,x1,x2,x3,New);
     inc(pos);
     end;

   for grupo:=1 to ng do for grupo1:=1 to ng do begin
     if grupo = grupo1 then begin
       SIGMABSXenon[canal,trozo,grupo] :=
         InterpolateTable(Q,NMAT,pos,a1,a2,a3,x1,x2,x3,New);
       SIGMABS[canal,trozo,grupo] := SIGMABS[canal,trozo,grupo] +
         DELTA*SIGMABSXenon[canal,trozo,grupo];
       end
     else
       SCATTERING[canal,trozo,grupo,grupo1] := SCATTERING[canal,trozo,grupo,grupo1] +
           DELTA*InterpolateTable(Q,NMAT,pos,a1,a2,a3,x1,x2,x3,New);
     inc(pos);
     end;

   for grupo:=1 to ng do begin
     NUSIGMAFIS[canal,trozo,grupo] := NUSIGMAFIS[canal,trozo,grupo] +
              DELTA*InterpolateTable(Q,NMAT,pos,a1,a2,a3,x1,x2,x3,New);
     inc(pos);
     end;

   for grupo:=1 to ng do begin
     SIGMAFIS[canal,trozo,grupo] := SIGMAFIS[canal,trozo,grupo] +
              DELTA*InterpolateTable(Q,NMAT,pos,a1,a2,a3,x1,x2,x3,New);
     inc(pos);
     end;

   end;  // agregar_Xenon

(* ---------------------------------------------------------------------- *)
 procedure agregar_Samario (N00: integer; DELTA: real);
 var
   grupo, grupo1, pos: integer;
 begin
   pos := N00+1;
   for grupo:=1 to ng do begin
     DIF[canal,trozo,grupo] := DIF[canal,trozo,grupo] +
       DELTA*InterpolateTable(Q,NMAT,pos,a1,a2,a3,x1,x2,x3,New);
     inc(pos);
     end;

   for grupo:=1 to ng do for grupo1:=1 to ng do begin
     if grupo = grupo1 then begin
       SIGMABSSamario[canal,trozo,grupo] :=
         InterpolateTable(Q,NMAT,pos,a1,a2,a3,x1,x2,x3,New);
       SIGMABS[canal,trozo,grupo] := SIGMABS[canal,trozo,grupo] +
         DELTA*SIGMABSSamario[canal,trozo,grupo];
       end
     else
       SCATTERING[canal,trozo,grupo,grupo1] := SCATTERING[canal,trozo,grupo,grupo1] +
           DELTA*InterpolateTable(Q,NMAT,pos,a1,a2,a3,x1,x2,x3,New);
     inc(pos);
     end;

   for grupo:=1 to ng do begin
     NUSIGMAFIS[canal,trozo,grupo] := NUSIGMAFIS[canal,trozo,grupo] +
              DELTA*InterpolateTable(Q,NMAT,pos,a1,a2,a3,x1,x2,x3,New);
     inc(pos);
     end;

   for grupo:=1 to ng do begin
     SIGMAFIS[canal,trozo,grupo] := SIGMAFIS[canal,trozo,grupo] +
              DELTA*InterpolateTable(Q,NMAT,pos,a1,a2,a3,x1,x2,x3,New);
     inc(pos);
     end;

   end;  // agregar_Samario

(* ----------------------------------------------------------------------- *)
 procedure imprimir_XS_por_trozo;
 var
   grupo, canal, trozo: integer;

 begin
  if not ArchivoSalidaAbierto then begin
    assign(F,ArchivoDeSalida);
    rewrite(F);
    ArchivoSalidaAbierto := true;
    end;

  writeln(F);
  writeln(F);
  for grupo:=1 to ng do begin
    writeln(F);
    writeln(F,'DIFUSION  GRUPO ',grupo);
    for canal:=1 to NCANALES do begin
      writeln(F);
      writeln(F,'canal ',canal);
      for trozo:=1 to NTROZOS do
        write(F,FloatToStrF(DIF[Canal,trozo,grupo],ffexponent,5,2),' ');
      end;
    end;

  writeln(F);
  writeln(F, 'ABSORCION');
  for grupo:=1 to ng do begin
     writeln(F);
     writeln(F,'ABSORCION  GRUPO ',grupo);
     for canal:=1 to NCANALES do begin
       writeln(F);
       writeln(F,'canal ',canal);
       for trozo:=1 to NTROZOS do
         write(F,FloatToStrF(SIGMABS[Canal,trozo,grupo],ffexponent,5,2),' ');
       end;
     end;

  writeln(F);
  writeln(F, 'NU SIGMAFIS');
  for grupo:=1 to ng do begin
      writeln(F);
      writeln(F,'NU SIGMAFIS  GRUPO ',grupo);
      for canal:=1 to NCANALES do begin
        writeln(F);
        writeln(F,'canal ',canal);
        for trozo:=1 to NTROZOS do
          write(F,FloatToStrF(NUSIGMAFIS[Canal,trozo,grupo],ffexponent,5,2),' ');
      end;
    end;

   Close(F);
   writeln ('Listo');
   Halt(0);

  end;   // imprimir_XS_por_trozo

(* --------------------------------------------------------------------- *)
  begin
    ttt0 := time();
    asignar_materiales;

 //  Con tablas de interpolación múltiple

    if ArchivoIntMultiple[1] <> '' then begin

      if not InterpolacionMultipleInicializada then begin
        for i:=1 to 6 do
          SetFileName(ArchivoIntMultiple[i],i,4);
        InterpolacionMultipleInicializada := true;
        end;

      for canal:=1 to NCANALES do
      for trozo:=1 to NTROZOS do begin
        P1 := BurnUp[canal,trozo];
        P2 := ConcXe[canal,trozo];
        P3 := TempRefr[canal,trozo];
        P4 := TComb[canal,trozo];
        P5 := DensRefr[canal,trozo];
        P6 := 0.0; //PPMBoro[canal,trozo]*1.0E-4;
        P7 := 0.0;
        Q := P1;

        fraccionado := NumMaterialBarra[canal,trozo] > 0;
        fraccion := FraccionBarra[canal,trozo];

        for caso := false to fraccionado do begin
          if caso then
            NMAT := NumMaterialBarra[canal,trozo]
          else
            NMAT := NumMaterial[canal,trozo];

          AllXS := Interpolate (P1,P2,P3,P4,P5,P6,P7,NMAT);

          pos := 1;
          for grupo:=1 to ng do begin
            if caso then
              DIF[canal,trozo, grupo] :=
                 (1.0-fraccion)*DIF[canal,trozo, grupo] + fraccion*ALLXS[pos]
            else
              DIF[canal,trozo, grupo] := ALLXS[pos];
            inc(pos);
            end;

          for grupo:=1 to ng do for grupo1:=1 to ng do begin
            if grupo = grupo1 then
              if caso then
                SIGMABS[canal,trozo,grupo] :=
                    (1.0-fraccion)*SIGMABS[canal,trozo,grupo] + fraccion*ALLXS[pos]
              else
                SIGMABS[canal,trozo,grupo] := ALLXS[pos]
            else
              if caso then
                SCATTERING[canal,trozo,grupo,grupo1] :=
                   (1.0-fraccion)*SCATTERING[canal,trozo,grupo,grupo1] + fraccion*ALLXS[pos]
              else
                SCATTERING[canal,trozo,grupo,grupo1] := ALLXS[pos];
            inc(pos);
            end;

          for grupo:=1 to ng do begin
            if caso then
              NUSIGMAFIS[canal,trozo, grupo] :=
                    (1.0-fraccion)*NUSIGMAFIS[canal,trozo, grupo] + fraccion*ALLXS[pos]
            else
              NUSIGMAFIS[canal,trozo, grupo] :=  ALLXS[pos];
            inc(pos);
            end;

          for grupo:=1 to ng do begin
            if caso then
              SIGMAFIS[canal,trozo, grupo] :=
                    (1.0-fraccion)*SIGMAFIS[canal,trozo, grupo] + fraccion*ALLXS[pos]
            else
              SIGMAFIS[canal,trozo, grupo] :=  ALLXS[pos];
            inc(pos);
            end;

          end;

        if ConXenon then begin
          table := AllTables[NMAT];
          N00 := 241;
          hallar_SIGMABS_Xenon (N00);
          end;

        end;

     end

 //  Con tablas de interpolación simple

    else
    for canal:=1 to NCANALES do
    for trozo:=1 to NTROZOS do begin
      fraccionado := NumMaterialBarra[canal,trozo] > 0;
      fraccion := FraccionBarra[canal,trozo];

      for caso := false to fraccionado do begin
        if caso then
          NMAT := NumMaterialBarra[canal,trozo]
        else
          NMAT := NumMaterial[canal,trozo];

        table := AllTables[NMAT];
        Q := BurnUp[canal,trozo];

        New:= true;
        pos := 2;
        for grupo:=1 to ng do begin
          val := InterpolateTable(Q,NMAT,pos,a1,a2,a3,x1,x2,x3,New);
          if val >= 0.0 then
            if caso then
              DIF[canal,trozo, grupo] := (1.0-fraccion)*DIF[canal,trozo, grupo] + fraccion/(3.0*val)
            else
              DIF[canal,trozo, grupo] := 1.0/(3.0*val);

          inc(pos);
          end;

        for grupo:=1 to ng do for grupo1:=1 to ng do begin
          val := InterpolateTable(Q,NMAT,pos,a1,a2,a3,x1,x2,x3,New);

          if caso then
            if grupo = grupo1 then
             SIGMABS[canal,trozo,grupo] :=
                 (1.0 - fraccion)*SIGMABS[canal,trozo,grupo] + fraccion*val
            else
             SCATTERING[canal,trozo,grupo,grupo1] :=
                (1.0 - fraccion)*SCATTERING[canal,trozo,grupo,grupo1] + fraccion*val
          else
            if grupo = grupo1 then
             SIGMABS[canal,trozo,grupo] := val
            else
             SCATTERING[canal,trozo,grupo,grupo1] := val;

          inc(pos);
          end;

        for grupo:=1 to ng do begin
          val := InterpolateTable(Q,NMAT,pos,a1,a2,a3,x1,x2,x3,New);
          if caso then
            NUSIGMAFIS[canal,trozo, grupo] :=
                (1.0 - fraccion)*NUSIGMAFIS[canal,trozo, grupo] + fraccion*val
          else
            NUSIGMAFIS[canal,trozo, grupo] := val;
          inc(pos);
          end;

        for grupo:=1 to ng do begin
          val := InterpolateTable(Q,NMAT,pos,a1,a2,a3,x1,x2,x3,New);
          if caso then
            SIGMAFIS[canal,trozo, grupo] :=
                (1.0 - fraccion)*SIGMAFIS[canal,trozo, grupo] + fraccion*val
          else
            SIGMAFIS[canal,trozo, grupo] := val;
          inc(pos);
          end;

        end;

      if ConReacoplamientoTermohidraulico then begin
        N00 := 41;
        DELTA := TComb[canal,trozo] - TComb0;
        agregar (N00, DELTA);

        N00 := 81;
        DELTA := TempRefr[canal,trozo] - TempRefr0;
        agregar (N00, DELTA);

        DELTA := DensRefr[canal,trozo] - 625.0;

        if DensRefr[canal,trozo] <= 625.0 then
          agregar (161,DELTA)
        else
          agregar (201,DELTA);
        end;

      if ConXenon then begin
        DELTA := ConcXe[canal,trozo] - ConcXe0;
        N00 := 241;
        agregar_Xenon (N00, DELTA);
        end;

      end;

 //  Con ambos tipos de tablas

    for canal:=1 to NCANALES do
    for trozo:=1 to NTROZOS do begin

      if ConSamario then begin
        DELTA := ConcSm[canal,trozo] - ConcSamario0;
        N00 := 281;
        agregar_Samario (N00, DELTA);
        end;

      if ConBoro then agregar_boro;
//        N00 := 161;
//        DELTA := PPMBoro[canal,trozo];
//        end;

      if (trozo = 6) or (trozo = 11) then agregar_separador;

      end;

    for canal:=1 to NCANALES do
    for trozo:=1 to NTROZOS do
    for grupo:=1 to ng do
      if DIF[canal,trozo, grupo] >= 0.0 then
        DIF[canal,trozo, grupo] := 1.0/(3.0*DIF[canal,trozo, grupo]);

  NUSIGMAFIS0 :=NUSIGMAFIS;

//  imprimir_XS_por_trozo;

 {   for canal:=1 to NTotalRegiones do begin
      j1 := reticula[canal,1];
      j2 := reticula[canal,2]-1;
      i1 := reticula[canal,3];
      i2 := reticula[canal,4]-1;

      for i:=i1 to i2 do
      for j:=j1 to j2 do
      CONTROL[1,i,j] := canal;
      end;
  }
    for iz:=NZMIN to NZMAX do begin
//    for iz:=1 to NZ do begin

      for canal:=1 to NCANALES do begin
        j1 := reticula[canal,1];
        j2 := reticula[canal,2]-1;
        i1 := reticula[canal,3];
        i2 := reticula[canal,4]-1;

        for i:=i1 to i2 do
        for j:=j1 to j2 do
        if iz < NREFL1 then begin
          pos := 1;
          for grupo:=1 to ng do begin
            D[iz,i,j,grupo] := XS_REFL_Inf[pos];
            inc(pos);
            end;

          for grupo:=1 to ng do for grupo1:=1 to ng do begin
            if grupo = grupo1 then
              SIGTOT[iz,i,j,grupo] := XS_REFL_Inf[pos]
            else
              SIGSCATT[iz,i,j,grupo,grupo1] :=  XS_REFL_Inf[pos];
            inc(pos);
            end;

          for grupo:=1 to ng do for grupo1:=1 to ng do
            SIGTOT[iz,i,j,grupo] := SIGTOT[iz,i,j,grupo]
                + SIGSCATT[iz,i,j,grupo,grupo1];

          for grupo:=1 to ng do begin
            NUSIGF[iz,i,j,grupo] :=  0.0;
            NUSIGF0[iz,i,j,grupo] :=  0.0;
            SIGF[iz,i,j,grupo] :=  0.0;
            end;

          end

        else if iz <= NREFL2 then begin
          trozo := (iz-NREFL1+1);
          D[iz,i,j] := DIF[canal,trozo];
          for grupo:=1 to ng do begin
            SIGTOT[iz,i,j,grupo] :=  SIGMABS[canal,trozo,grupo];
            for grupo1:=1 to ng do
              SIGTOT[iz,i,j,grupo] := SIGTOT[iz,i,j,grupo]
                + SCATTERING[canal,trozo,grupo,grupo1];

            end;

          SIGSCATT[iz,i,j] :=  SCATTERING[canal,trozo];

          NUSIGF[iz,i,j] :=  NUSIGMAFIS[canal,trozo];
          NUSIGF0[iz,i,j] :=  NUSIGMAFIS0[canal,trozo];
          SIGF[iz,i,j] :=  SIGMAFIS[canal,trozo];
          end

        else begin
          pos := 1;
          for grupo:=1 to ng do begin
            D[iz,i,j,grupo] := XS_REFL_Sup[pos];
            inc(pos);
            end;
          for grupo:=1 to ng do for grupo1:=1 to ng do begin
            if grupo = grupo1 then
              SIGTOT[iz,i,j,grupo] := XS_REFL_Sup[pos]
            else
              SIGSCATT[iz,i,j,grupo,grupo1] :=  XS_REFL_Sup[pos];
            inc(pos);
            end;

          for grupo:=1 to ng do for grupo1:=1 to ng do
            SIGTOT[iz,i,j,grupo] := SIGTOT[iz,i,j,grupo]
                + SIGSCATT[iz,i,j,grupo,grupo1];

          for grupo:=1 to ng do begin
            NUSIGF[iz,i,j,grupo] :=  0.0;
            NUSIGF0[iz,i,j,grupo] :=  0.0;
            SIGF[iz,i,j,grupo] :=  0.0;
            end;

          end

        end;

      end;

    for iz:=NZMIN to NZMAX do begin
//    for iz:=1 to NZ do begin

      for canal:=62 to NTotalRegiones do begin
        j1 := reticula[canal,1];
        j2 := reticula[canal,2]-1;
        i1 := reticula[canal,3];
        i2 := reticula[canal,4]-1;

        if canal <= NRegionesReactor then
        for i:=i1 to i2 do
        for j:=j1 to j2 do begin

          if iz < NREFL1 then begin
            pos := 1;
            for grupo:=1 to ng do begin
              D[iz,i,j,grupo] := XS_REFL_Rad[pos];
              inc(pos);
              end;
            for grupo:=1 to ng do for grupo1:=1 to ng do begin
              if grupo = grupo1 then
                SIGTOT[iz,i,j,grupo] := XS_REFL_Rad[pos]
              else
                SIGSCATT[iz,i,j,grupo,grupo1] :=  XS_REFL_Rad[pos];
              inc(pos);
              end;

            for grupo:=1 to ng do for grupo1:=1 to ng do
              SIGTOT[iz,i,j,grupo] := SIGTOT[iz,i,j,grupo]
                  + SIGSCATT[iz,i,j,grupo,grupo1];

            end

          else if iz > NREFL2 then begin
            pos := 1;
            for grupo:=1 to ng do begin
              D[iz,i,j,grupo] := XS_REFL_Rad[pos];
              inc(pos);
              end;
            for grupo:=1 to ng do for grupo1:=1 to ng do begin
              if grupo = grupo1 then
                SIGTOT[iz,i,j,grupo] := XS_REFL_Rad[pos]
              else
                SIGSCATT[iz,i,j,grupo,grupo1] :=  XS_REFL_Rad[pos];
              inc(pos);
              end;

            for grupo:=1 to ng do for grupo1:=1 to ng do
              SIGTOT[iz,i,j,grupo] := SIGTOT[iz,i,j,grupo]
                  + SIGSCATT[iz,i,j,grupo,grupo1];

            end

          else begin
            pos := 1;
            for grupo:=1 to ng do begin
              D[iz,i,j,grupo] := XS_REFL_Rad[pos];
              inc(pos);
              end;
            for grupo:=1 to ng do for grupo1:=1 to ng do begin
              if grupo = grupo1 then
                SIGTOT[iz,i,j,grupo] := XS_REFL_Rad[pos]
              else
                SIGSCATT[iz,i,j,grupo,grupo1] :=  XS_REFL_Rad[pos];
              inc(pos);
              end;

            for grupo:=1 to ng do for grupo1:=1 to ng do
              SIGTOT[iz,i,j,grupo] := SIGTOT[iz,i,j,grupo]
                  + SIGSCATT[iz,i,j,grupo,grupo1];
            end;

            for grupo:=1 to ng do begin
              NUSIGF[iz,i,j,grupo] :=  0.0;
              NUSIGF0[iz,i,j,grupo] :=  0.0;
              SIGF[iz,i,j,grupo] :=  0.0;
            end;
          end

        else
        for i:=i1 to i2 do
        for j:=j1 to j2 do begin
          for grupo:=1 to ng do begin
            D[iz,i,j,grupo] := 0.0;
            SIGTOT[iz,i,j,grupo] :=  0.0;
            for grupo1:=1 to ng do
              SIGSCATT[iz,i,j,grupo,grupo1] :=  0.0;
            NUSIGF[iz,i,j,grupo] :=  0.0;
            SIGF[iz,i,j,grupo] :=  0.0;
            end
          end

        end;

      end;

//    ttt0 := time-ttt0;
//    writeln (F);
//    write (F,'  Tiempo TOTAL (CPU) para las XS: ',ttt0*86400.0:10:3,' seg');
    end;

//********************************************************************
//*  Lectura de todos los datos de entrada de la unidad definida     */
//*  en ArchivoDeEntrada. Previamente se eliminan los registros en   */
//*  blanco y los comentarios y se pasa todo a caracteres mayúsculos */
//********************************************************************/

 procedure ReadDATA;
 var
   canal,trozo,i,j,nt: integer;
   S: string;
   Entrada: text;
   ArchivoTablas: Text;
   Intercalado, salir: boolean;

 (* --------------------------------------- *)
 function E(ss: string): boolean;
 begin
   Result := pos(ss,S) > 0;
   end;

 (* --------------------------------------- *)
 function trim(ss: string) : string;
 var
   n: integer;

 begin

   repeat
     n := pos(' ',ss);
     if n > 0 then delete(ss,n,1);
     until n = 0;

   Result := ss;
   end;

 (* --------------------------------------- *)
 procedure Uppercase(var SS: string);
 var
   i: integer;
   c: char;

 begin
   for i:=1 to length(SS) do begin
     c := SS[i];
     if (c in ['a'..'z']) then c := Upcase(c)
     else if (c in ['á','Á']) then c := 'A'
     else if (c in ['é','É']) then c := 'E'
     else if (c in ['í','Í']) then c := 'I'
     else if (c in ['ó','Ó']) then c := 'O'
     else if (c in ['ú','Ú']) then c := 'U';

     SS[i] := c;
     end;

   end;

 (* --------------------------------------- *)
 begin
   Intercalado := false;
   Salir := false;

   AssignFile(Entrada, ArchivoDeEntrada);
   Reset(Entrada);

   if ArchivoDeSalida <> '' then begin
     if ArchivoSalidaAbierto then
       Close(F);

     AssignFile(F,ArchivoDeSalida);
     Rewrite(F);
     ArchivoSalidaAbierto := true;

     writeln(F);
     writeln(F,'*********************************************************************');
     writeln(F);
     writeln(F,'   Se leyeron los datos del archivo: ',ArchivoDeEntrada);
     writeln(F,'   Se graba la salida en el archivo: ',ArchivoDeSalida);
     writeln(F);
     writeln(F,'*********************************************************************');
     writeln(F);
     end;

   while (not eof(Entrada)) and (not salir) do begin
     readln(Entrada,S);

     uppercase(S);

     if (S = '') or (S[1] = '!') then

     else if S[1] = '*' then
       if ArchivoDeSalida <> '' then
         writeln(F,S)
       else

     else if E('FIN') and E('DATOS') then
       salir := true

     else if E('TFINAL') then
       readln (Entrada,TFINAL)

     else if E('OPCIONES') then begin
       ConXenon := E('XENON');
       ConSamario := E('SAMARIO');
       ConReacoplamientoTermohidraulico := E('TERMOHIDRAULIC') and E('REACOPL');
       ConBarras := E('BARRAS');
       ConQuemado := E('QUEMADO');
       cinetica := E('CINETICA');
       variacion_lineal := E('LINEAL');
       if cinetica then MetodoAdiabatico := not E('DIRECTO');
       ConBoro := E('BORO');
//       ciclo := E('CICLO');
       end

     else if E('GUARDAR') and E('ESTADO') then begin
        readln (Entrada,S);
        SaveFileName := trim(S);
	GuardarEstado := true;
        end

    else if E('LEER') and E('ESTADO') then begin
      readln(Entrada,S);
      RetrieveFileName := trim(S);
      LeerEstado := true;
      end

     else if E('MATRIZ') and E('RESPUESTA') then ConMatrizDeRespuesta := true

     else if E('MAXIMO') and E('ITERACIONES') then
        readln(Entrada,MaxIteraciones)

     else if E('IMPRIMIR') then begin
       ImprimirIteraciones := E('ITERACIONES');
       ImprimirNumIteraciones := E('NUMERO');
       ImprimirDistribucion := E('DISTRIBUCION');
       ImprimirTempDens := E('TEMPERATURA');
       ImprimirPorCanal := E('POR') and E('CANAL');
       end

     else if E('PUMA') and E('FACTOR') then begin
       ParaCompararConPUMA := true;
       readln (Entrada,FactPuma);
       end

     else if E('COEFICIENTE') and E('SOBRERRELAJACION') then
       readln (Entrada,BETA0)

     else if E('PASO') and E('TIEMPO') then
       readln (Entrada,DELTAT)

     else if E('PASO') and E('IMPRESION') then
       readln (Entrada,PasoDeImpresion)

     else if E('PRECISION') then begin

       if E('CINETICO') then
         readln (Entrada, PrecisionCinetico)
       else
	 readln (Entrada,PrecisionEstacionario)

       end

     else if E('INTERPOLACION') and E('LINEAL') then
        interpolacion_lineal := true

     else if E('INTERPOLACION') and E('MULTIPLE') then
        for i:=1 to 6 do
        readln(Entrada,ArchivoIntMultiple[i])

     else if E('POTENCIA') and E('FISION') then
       readln (Entrada,PPorFision)

     else if E('INICIAL') and E('VALORES') and E('CONSTANTES') then begin
       readln (Entrada,ConcXE0, TComb0, DensRefr0, TempRefr0);
       for canal:=1 to NCANALES do
       for trozo:=1 to NTROZOS do begin
         DensRefr[canal,trozo] := DensRefr0;
         TComb[canal,trozo] := Tcomb0;
         ConcXe[canal,trozo] := ConcXe0;
         TempRefr[canal,trozo] := TempRefr0;
         end;

       end

     else if E('TABLAS') then begin
       readln(Entrada,S);
       Assign(ArchivoTablas,S);

       Reset(ArchivoTablas);

       for nt:=1 to NMaxTablas do begin
         Readln(ArchivoTablas,NumFilas[nt]);

         for i:=1 to NumFilas[nt] do
         for j:=1 to LongTableLine0 do
           read (ArchivoTablas, AllTables[nt,i,j]) ;
         end;

       Close (ArchivoTablas);
       end

     else if E('REFLECTOR') then begin
       for j:=1 to 40 do read (Entrada, XS_REFL_Inf[j]) ;
       for j:=1 to 40 do read (Entrada, XS_REFL_Sup[j]) ;
       for j:=1 to 40 do read (Entrada, XS_REFL_Rad[j]) ;
       end

     else if E('INTERCALADO') then
       intercalado := true

     else if E('MATERIALES') then begin
         for canal:=1 to NCANALES do
         for trozo:=1 to NTROZOS do
           read (Entrada,NumMaterial[canal,trozo]);
       end

     else if E('TEMPERATURA') and E('COMBUSTIBLE') then begin
       for canal:=1 to NCANALES do begin
         if intercalado then
           readln(Entrada);

         for trozo:=1 to NTROZOS do
           read (Entrada,TCOMB[RENUM_INV[canal],trozo]);
         readln(Entrada);
         end;

       end

     else if E('TEMPERATURA') and E('REFRIGERANTE') then begin
       for canal:=1 to NCANALES do begin
         if intercalado then
           readln(Entrada);

         for trozo:=1 to NTROZOS do
           read (Entrada,TempRefr[RENUM_INV[canal],trozo]);
         readln(Entrada);
         end;
       end

     else if E('DENSIDAD') and E('REFRIGERANTE') then begin
       for canal:=1 to NCANALES do begin
         if intercalado then
           readln(Entrada);

         for trozo:=1 to NTROZOS do
           read (Entrada,DensRefr[RENUM_INV[canal],trozo]);
         readln(Entrada);
         end;
       end

     else if E('QUEMADO') then
       for canal:=1 to NCANALES do begin
         if intercalado then
           readln(Entrada);

         for trozo:=1 to NTROZOS do begin
           read (Entrada,BurnUp[RENUM_INV[canal],trozo]);
           end;
         readln(Entrada);
         end;

   end;  //  while

 if not ImprimirIteraciones and not ImprimirNumIteraciones
   and not ImprimirDistribucion and not ImprimirPorCanal then begin
    close(F);
    ArchivoSalidaAbierto := false;
    end;

 CloseFile (Entrada);
 end;  // ReadDATA

(*******************************************************************************)
 procedure inicializar;
 var
   canal,trozo, prec,grupo, i, j, iz, i1, i2, j1, j2: integer;

 begin
//   ArchivoSalidaAbierto := false;
  ConSamario := false;
  ConXenon := false;                   ConReacoplamientoTermohidraulico := false;
  ImprimirIteraciones := false;        ImprimirNumIteraciones := false;
  ConBarras := false;                  ConQuemado := false;
  ParaCompararConPUMA := false;        ImprimirIteraciones := false;
  ImprimirNumIteraciones := false;     ImprimirDistribucion := false;
  ImprimirPorCanal := false;
  UNIDAD := 1.0;                       FactPuma := 1.0;
  PrecisionEstacionario := 1.0E-6;     PrecisionCinetico := 1.0E-6;
  PuedeImprimir := true;               PPorFision := 1.0;
  MaxIteraciones := 600;               variacion_lineal := false;
  InterpolacionMultipleInicializada := false;
  ConMatrizDeRespuesta := false;
  BETA0 := 1.6;                        interpolacion_lineal := false;

  ConcSamario0 := 0.0; //4.4004E+16; // 1.0856E16;
  ConcXE0 := 1.3591E+15;
  TComb0 := 455.0;
  DensRefr0 := 625.0;
  TempRefr0 := 305.0;

  for i:=1 to 6 do ArchivoIntMultiple[i] := '';

  PasoDeImpresion := 0.0;
  ProximaImpresion := 0.0;
  cinetica := true;
  ciclo := false;
  GuardarEstado := false;
  LeerEstado := false;
  DeltaT := 0.1;
  TFINAL := 0.0;
  NumIteracionesEstacionarias := 15;
  TTime := 0.0;

  BetaNuclearTotal := 0.0;
  for i:=1 to nprec do
    BetaNuclearTotal := BetaNuclearTotal + BetaNuclear[i];

   for grupo:=1 to ng do begin
     espectro_estatico[grupo] := (1.0 - BetaNuclearTotal)*espectro_de_fision[grupo];

     for prec := 1 to nprec do
       espectro_estatico[grupo] := espectro_estatico[grupo] + chis[prec,grupo]*BetaNuclear[prec];

     end;

  for canal:=1 to NCANALES do
  for trozo:=1 to NTROZOS do begin
    DensRefr[canal,trozo] := DensRefr0;
    TComb[canal,trozo] := Tcomb0;
    ConcXe[canal,trozo] := ConcXe0;
    ConcSm[canal,trozo] := ConcSamario0;
    TempRefr[canal,trozo] := TempRefr0;
    PPMBoro[canal,trozo] := 0.0;
    BurnUp[canal,trozo] := 0.0;
    end;

  for canal:=1 to NCANALES do
    for trozo:=1 to NTROZOS do
      NumMaterial[canal,trozo] := 1;

  XS_REFL_Inf := XS_REFL_Inf0;
  XS_REFL_Sup := XS_REFL_Sup0;
  XS_REFL_Rad := XS_REFL_Rad0;

  ReadDATA;

  NZMIN := 1;
  NZMAX := NZ;

  FracPotResidual := 0.0;
  for i:=1 to NCOF do
    FracPotResidual := FracPotResidual + cof[i]/cofexp[i];

  KE0 := 1.0;
  KE_ANT := 1.0;

  for canal:=1 to NCANALES do
  for trozo:=1 to NTROZOS do
  for grupo:=1 to ng do
    flujo[canal,trozo,grupo] := 1.0;

  for iz:=1 to NZ do
  for i:=1 to NY do
  for j:=1 to NX do
  for grupo:=1 to ng do begin
    FI[iz,i,j,grupo] := 1.0;
    FIA[iz,i,j,grupo] := 1.0;
    end;

  VOLCeldilla := DELTAZ*ALTURA*ALTURA/sqrt(3);
  VolTrozo := 6.0*VOLCeldilla;

  for i:=1 to 6 do
    BetaLambda[i] := BetaNuclear[i]/lambda[i];

  end;

 //**************************************************************************/
//*  Se calculan las potencias residuales para cada canal y trozo          */
//*  para el estado estacionario (DT = 0) y para el no estacionario        */
//**************************************************************************/

 procedure potencias_residuales(DT: real);
 var
   sp1: CoreMatrix;
   k,canal,trozo: integer;
   x,y,xy: real;

 begin
   if DT = 0.0 then begin

     for canal:=1 to NCANALES do
     for trozo:=1 to NTROZOS do begin
       PotRes[canal,trozo] := PotInst[canal,trozo]*FracPotResidual;
       for k:=1 to NCOF do PAux[k,canal,trozo] := 0.0;
       PotInstAnt[canal,trozo] := PotInst[canal,trozo];
       end;

     end

   else begin
     for canal:=1 to NCANALES do
     for trozo:=1 to NTROZOS do
       sp1[canal,trozo] := PotInst[canal,trozo] - PotInstAnt[canal,trozo];

     for k:=1 to ncof do begin
       x := DT*cofexp[k];

       if abs(x) < 0.1 then
         xy := x*(1+x*(-0.5+x*(0.166667+x*(-0.0416667+x*0.00833333))))/cofexp[k]
       else
         xy := (1-exp(-x))/cofexp[k];

       y := exp(-x);

       for canal:=1 to NCANALES do
       for trozo:=1 to NTROZOS do begin
         paux[k,canal,trozo] := paux[k,canal,trozo] + cof[k]*sp1[canal,trozo];
         PotRes[canal,trozo] := PotRes[canal,trozo] + paux[k,canal,trozo]*xy;
         paux[k,canal,trozo] := paux[k,canal,trozo]*y;
         end;

       end;

     for canal:=1 to NCANALES do
     for trozo:=1 to NTROZOS do
       PotInstAnt[canal,trozo] := PotInst[canal,trozo];

     end;

   PotResidual := 0.0;
   for canal:=1 to NCANALES do
   for trozo:=1 to NTROZOS do begin
     PotResidual := PotResidual + PotRes[canal,trozo];
     POT[canal,trozo] := PotRes[canal,trozo] + PotInst[canal,trozo];
     end;

   end;   (* potencias_residuales  *)

//**************************************************************************/
//*  Calcula la distribución de flujos por canales y trozos a partir de    */
//*  la distribución de flujos para celdillas                              */
//**************************************************************************/

 procedure PasarFlujos;
 var
   iz,grupo,i,j,i1,i2,j1,j2: integer;
   canal,trozo: integer;

 begin

   for grupo:=1 to ng do
   for canal:=1 to NCANALES do
   for trozo:=1 to NTROZOS do begin
     flujo[canal,trozo,grupo] := 0.0;
     end;

   for iz:=NREFL1 to NREFL2 do begin
     trozo := (iz-NREFL1+1);
//     if NZ > 20 then
//       if NTROZOS = 10 then
//         trozo := (iz-NREFL1+2) div 2
//       else
//         trozo := (iz-NREFL1+1)
//     else
//       trozo := (iz-NREFL1+1);

     for canal:=1 to NCANALES do begin
       j1 := reticula[canal,1];
       j2 := reticula[canal,2]-1;
       i1 := reticula[canal,3];
       i2 := reticula[canal,4]-1;

       for i:=i1 to i2 do
       for j:=j1 to j2 do
       for grupo:=1 to ng do
         flujo[canal,trozo,grupo] := flujo[canal,trozo,grupo]
              + FI[iz,i,j,grupo]*VolCeldilla;

       end;

     end;

   for canal:=1 to NCANALES do
   for trozo:=1 to NTROZOS do
   for grupo:=1 to ng do
     flujo[canal,trozo,grupo] := flujo[canal,trozo,grupo]/VOLTROZO;

   end;

//**************************************************************************/
//*  Calcula la distribución de potencias para cada canal y trozo,         */
//*  potencias totales por canal, potencias y factores de forma por sector,*/
//*  factores de asimetría, factor de forma de la potencia                 */
//**************************************************************************/

 procedure potencias(POWER,DT: real);
 var
   canal,trozo,grupo,sxt,iz,i,j: integer;
   ACC, FNORM, EX1, EX2: real;
   NCSexto: array[1..6] of integer;
   ACCSexto: array[1..6,1..2] of real;

 begin
  PasarFlujos;

  ACC := 0.0;

  for canal:=1 to NCANALES do
  for trozo:=1 to NTROZOS do begin
    PotInst[canal,trozo] := 0.0;

    for grupo:=1 to ng do
      PotInst[canal,trozo] := PotInst[canal,trozo] +
       flujo[canal,trozo,grupo]*SIGMAFIS[canal,trozo,grupo]
          *VOLTROZO*PPorFision*1.6E-19;

    ACC := ACC + PotInst[canal,trozo];
    end;

  if DT = 0.0 then FNORM := POWER/ACC/(1.0 + FracPotResidual)
              else FNORM := 1.0;

  for canal:=1 to NCANALES do
  for trozo:=1 to NTROZOS do begin
    for grupo:=1 to ng do
      flujo[canal,trozo,grupo] := flujo[canal,trozo,grupo]*FNORM;
    PotInst[canal,trozo] := PotInst[canal,trozo]*FNORM;
    end;

  for iz:=1 to NZ do
  for i:=1 to NY do
  for j:=1 to NX do
  for grupo:=1 to ng do
    FI[iz,i,j,grupo] := FI[iz,i,j,grupo]*FNORM;

  FIA := FI;

  potencias_residuales(DT);

  PMax := 0.0;
  PotMax := 0.0;
  PotTotal := 0.0;
  CanalPMax := 0;
  TrozoPMax := 0;
  CanalPotMax := 0;

  for canal:=1 to NCANALES do begin
    ACC := 0.0;

    for trozo:=1 to NTROZOS do begin

      if POT[canal,trozo] > PMax then begin
        PMax := POT[canal,trozo];
        CanalPMax := canal;
        trozoPMax := trozo;
        end;

      ACC := ACC + POT[canal,trozo];

      end;

    potcan[canal] := ACC;
    PotTotal := PotTotal + ACC;

    if potcan[canal] > PotMax then begin
      PotMax := potcan[canal];
      CanalPotMax := canal;
      end;

    end;

  FForma := PMax/(PotTotal/(NCANALES*NTROZOS));

  end;

//**************************************************************************/
//*  Calcula la distribución de Xenón para régimen estacionario (DT = 0)   */
//*  y no estacionario, además de su valor medio                           */
//**************************************************************************/

 procedure XENON (DT: real);
 const
    ani = 2.93e-5;
    anx = 2.10e-5;

 var
   canal,trozo,grupo: integer;
   ProdIod, ProdXe, AbsXe, a, xn, xa, xan: real;
   concxe_ant, conci_ant: real;
   TFISS, TABS: real;

 const
   YieldXE = 0.002576 ;
   YieldIodo  = 0.062810 ;

//   SIod1: real = 1.2645E-04;
//   SIod2: real = 2.5753E-03;
//   SXe1: real = 3.5894E-06;
//   SXe2: real = 1.0559E-04;
//   SigAXe1: real = 1.3219E+02;
//   SigAXe2: real = 1.5072E+06;

 begin
 (*   estacionario  *)
   ConcXeMedia := 0.0;

   if DT = 0.0 then

     for canal:=1 to NCANALES do
     for trozo:=1 to NTROZOS do begin

       TFISS := 0.0;
       for grupo:=1 to ng do
         TFISS := TFISS + flujo[canal,trozo,grupo]* SIGMAFIS[canal,trozo, grupo];

       ProdIod := TFISS*YieldIodo;
       ProdXe := TFISS*YieldXE;

       AbsXe := 0.0;
       for grupo:=1 to ng do
         AbsXe := AbsXe + flujo[canal,trozo,grupo]* SIGMABSXenon[canal,trozo, grupo];

       ConcI[canal,trozo] := ProdIod/ani;
       ConcXE[canal,trozo] := (ProdIod + ProdXe)/(AbsXe+anx) ;
       ConcXeMedia := ConcXeMedia + ConcXE[canal,trozo];
       end

 (*  no estacionario  *)

   else
     for canal:=1 to NCANALES do
     for trozo:=1 to NTROZOS do begin

       TFISS := 0.0;
       for grupo:=1 to ng do
         TFISS := TFISS + flujo[canal,trozo,grupo]* SIGMAFIS[canal,trozo, grupo];

       ProdIod := TFISS*YieldIodo;
       ProdXe := TFISS*YieldXE;

       AbsXe := 0.0;
       for grupo:=1 to ng do
         AbsXe := AbsXe + flujo[canal,trozo,grupo]* SIGMABSXenon[canal,trozo, grupo];

       a := AbsXe+anx;
       xn := exp(-ani*DT);
       xa := exp(-a*DT);
       xan := (xn-xa)/(a-ani) ;
       concxe_ant := ConcXe[canal,trozo];
       conci_ant := ConcI[canal,trozo];

       concxe[canal,trozo]:= concxe_ant*xa + conci_ant*ani*xan
             + (ProdIod+ProdXe)*(1-xa)/a-ProdIod*xan;
       conci[canal,trozo] := conci_ant*xn + (1-xn)*ProdIod/ani;
       ConcXeMedia := ConcXeMedia + ConcXE[canal,trozo];
       end;

    ConcXeMedia := ConcXeMedia/(NCANALES*NTROZOS);
    if DT = 0.0 then ConcXenon0 := ConcXeMedia;
    end;

//**************************************************************************/
//*  Calcula la distribución de Xenón para régimen estacionario (DT = 0)   */
//*  y no estacionario, además de su valor medio                           */
//**************************************************************************/

 procedure SAMARIO (DT: real);
 const
    anp = 3.625E-6;
    ans = 0.0;

 var
   canal,trozo,grupo: integer;
   ProdPt, ProdSm, AbsSm, a, xn, xa, xan: real;
   concsm_ant, concpt_ant: real;
   TFISS, TABS: real;

 const
   YieldPt = 0.01081;
   YieldSm = 5.0420E-14;

 begin
 (*   estacionario  *)
   ConcSmMedia := 0.0;

   if DT = 0.0 then

     for canal:=1 to NCANALES do
     for trozo:=1 to NTROZOS do begin
       TFISS := 0.0;
       for grupo:=1 to ng do
         TFISS := TFISS + flujo[canal,trozo,grupo]* SIGMAFIS[canal,trozo, grupo];

       ProdPt := TFISS*YieldPt;

       AbsSm := 0.0;
       for grupo:=1 to ng do
         AbsSm := AbsSm + flujo[canal,trozo,grupo]* SIGMABSSamario[canal,trozo, grupo];

       ConcPt[canal,trozo] := ProdPt/anp;
       ConcSm[canal,trozo] := ProdPt/AbsSm ;
       ConcSmMedia := ConcSmMedia + ConcSm[canal,trozo];
//       ConcSm[canal,trozo] := ConcSm[canal,trozo] + ConcPt[canal,trozo];
       end

 (*  no estacionario  *)

   else
     for canal:=1 to NCANALES do
     for trozo:=1 to NTROZOS do begin

       TFISS := 0.0;
       for grupo:=1 to ng do
         TFISS := TFISS + flujo[canal,trozo,grupo]* SIGMAFIS[canal,trozo, grupo];

       ProdPt := TFISS*YieldPt;
       ProdSm := 0.0;

       AbsSm := 0.0;
       for grupo:=1 to ng do
         AbsSm := AbsSm + flujo[canal,trozo,grupo]* SIGMABSSamario[canal,trozo, grupo];

       a := AbsSm+ans;
       xn := exp(-anp*DT);
       xa := exp(-a*DT);
       xan := (xn-xa)/(a-anp) ;
       concsm_ant := ConcSm[canal,trozo];
       concPt_ant := ConcPt[canal,trozo];

       concsm[canal,trozo]:= concsm_ant*xa + concpt_ant*anp*xan
             + (ProdPt+ProdSm)*(1-xa)/a-ProdPt*xan;
       concpt[canal,trozo] := concpt_ant*xn + (1-xn)*Prodpt/anp;
       ConcSmMedia := ConcSmMedia + ConcSm[canal,trozo];
       end;

    ConcSmMedia := ConcSmMedia/(NCANALES*NTROZOS);
    if DT = 0.0 then ConcSamario0 := ConcSmMedia;
    end;

//*******************************************************/
//*  Imprime en las salidas pedidas                     */
//*******************************************************/
 procedure salidas (TTime, DT: real);
 var
   canal,trozo: integer;
   ACC: real;

(* ----------------------------------------------------- *)
 begin
  if not ArchivoSalidaAbierto then begin
    assign(F,ArchivoDeSalida);
    rewrite(F);
    ArchivoSalidaAbierto := true;
    end;

  if (not PuedeImprimir) then exit;

  writeln(F);
  writeln(F);
  writeln(F,'*****************************************************************');
  writeln (F,'Instante = ',TTime/UNIDAD:1:4);
  writeln(F);
  writeln (F,'     Potencia Total = ',PotTotal:1:4,'    ');
  writeln (F,'Potencia neutrónica = ',PotTotal-PotResidual:1:4);
  writeln (F,'  Potencia Residual = ',PotResidual:1:4);
  writeln (F,'    Factor de Forma = ',FForma:1:3);

  if (DT > 0.0) and cinetica then begin
    writeln(F);
    writeln (F,'Derivada Logarítmica = ',
                    100.0*DerivadaLogaritmica:10:3,'%');
    writeln(F);
    end;

  if ConXenon then begin
    writeln (F,'   Conc Xenon Media = ',FloatToStrF(ConcXEMedia,ffexponent,5,2));
    writeln (F,'  Conc Xen Relativa = ',ConcXEMedia/ConcXenon0:1:5);
    end;

  if ConSamario then begin
    writeln (F,'   Conc Samario Media = ',FloatToStrF(ConcSmMedia,ffexponent,5,2));
    writeln (F,'  Conc Samario Relativa = ',ConcSmMedia/ConcSamario0:1:5);
    end;

  writeln (F,'Pot Esp Máxima: ',
      PotTotal/(NCANALES*NTROZOS)*FForma/VOLTROZO*1.0E6:1:3,
          '   en el canal: ',RENUM[canalPMax]:3,'  trozo: ',TrozoPMax:3);

  if DT = 0.0 then begin
    writeln (F,'     KE (estático) = ',KE0:1:5);
    writeln (F,'     RO (estático) = ',((KE0-1.0)/KE0*100000.0):1:1,' PCM');
    writeln (F,'  React (estática) = ',RO*100000.0:1:1,' PCM');
    end

  else begin
    writeln (F,'       Reactividad = ',RO*100000.0:1:1,' PCM');
    writeln (F,'        RO/BETAEFF = ',RO/BetaNuclearTotal:1:3,' $');
    end;

//  writeln (F);
//  writeln (F,'Posiciones de los bancos de barras:');

  if ImprimirPorCanal then begin
    writeln (F);
    writeln (F);
    writeln(F,'****************************************************');
    writeln(F,'POTENCIAS POR CANAL');
    writeln (F);

    for canal:=1 to NCANALES do begin
      write (F,'CAN ',NomCan[canal],' ',potcan[RENUM_INV[canal]]:1:3,'  ');
      if canal mod 6 = 0 then writeln(F);
      end;

   end;

 if ImprimirDistribucion then begin
   writeln(F);
   writeln(F);
   writeln(F,'****************************************************');
   writeln(F,'DISTRIBUCIÓN ESPACIAL DE POTENCIA');
   writeln(F);
   for canal:=1 to NCANALES do begin
     writeln(F);
     writeln (F,'canal ',NomCan[canal],'  ',potcan[RENUM_INV[canal]]:1:2);
     for trozo:=1 to NTROZOS do begin
        write(F,POT[RENUM_INV[canal],trozo]/VOLTROZO*1.0E6:8:2,' ');
        if trozo mod 10 = 0 then writeln(F);
        end;
     end;

   ACC := 0.0;

   writeln(F);
   writeln(F);
   writeln(F,'****************************************************');
   writeln(F,'DISTRIBUCIÓN DE QUEMADO');

   for canal:=1 to NCANALES do begin
     writeln(F);
     writeln (F,'canal ',NomCan[canal]);
     for trozo:=1 to NTROZOS do begin
        write(F,BurnUp[RENUM_INV[canal],trozo]:8:1,' ');
        ACC := ACC + BurnUp[RENUM_INV[canal],trozo];
        if trozo mod 10 = 0 then writeln(F);
        end;
     end;

   writeln(F);
   writeln(F);
   write(F,'QUEMADO MEDIO :',ACC/(NCANALES*NTROZOS):8:1);
   writeln(F);

   if ImprimirTempDens then begin
     writeln(F);
     writeln(F);
     writeln(F,'****************************************************');
     writeln(F,'DISTRIBUCIÓN DE TEMPERATURA DEL COMBUSTIBLE');

     for canal:=1 to NCANALES do begin
       writeln(F);
       writeln (F,'canal ',NomCan[canal]);
       for trozo:=1 to NTROZOS do begin
          write(F,TComb[RENUM_INV[canal],trozo]:8:2,' ');
          if trozo mod 10 = 0 then writeln(F);
          end;
       end;

     writeln(F);
     writeln(F);
     writeln(F,'****************************************************');
     writeln(F,'DISTRIBUCIÓN DE TEMPERATURA DEL RFRIGERANTE');

     for canal:=1 to NCANALES do begin
       writeln(F);
       writeln (F,'canal ',NomCan[canal]);
       for trozo:=1 to NTROZOS do begin
          write(F,TempRefr[RENUM_INV[canal],trozo]:8:2,' ');
          if trozo mod 10 = 0 then writeln(F);
          end;
       end;

     writeln(F);
     writeln(F);
     writeln(F,'****************************************************');
     writeln(F,'DISTRIBUCIÓN DE DENSIDAD DEL REFRIGERANTE');

     for canal:=1 to NCANALES do begin
       writeln(F);
       writeln (F,'canal ',NomCan[canal]);
       for trozo:=1 to NTROZOS do begin
          write(F,DensRefr[RENUM_INV[canal],trozo]:8:2,' ');
          if trozo mod 10 = 0 then writeln(F);
          end;
       end;

     end;

   writeln(F);
   writeln(F);
   writeln(F,'****************************************************');
   writeln(F,'DISTRIBUCIÓN DE MATERIALES');
   writeln(F);

   for canal:=1 to NCANALES do begin
     writeln (F,'canal ',NomCan[canal]);
     for trozo:=1 to NTROZOS do begin
        write(F,NumMaterial[RENUM_INV[canal],trozo]:4,' ');
        end;
      writeln(F);
     end;

    end;

  writeln(F);
  end;

//**************************************************************************/
//*  Resuelve el problema de autovalores de flujo y K efectivo para la     */
//*  matriz hallada en estado estacionario                                 */
//**************************************************************************/

 procedure Iteraciones (TipoCalculo: TipoCalculoType; precision: real);
 var
   i,j,iz, iteracion, NNN, n, NP, MaxIter, grupo, grupo1: integer;
   EPS0, EX, EPSKE: real;
   FLU : array[1..ng] of real;
   Diagonal, TermIzquierdo: real;
   OMEGA, EPSOMEGA, OMEGA0, OMEGAX: real;
   KE, EPS, ALFA, ACC, ACC0: real;
   TTT00: TDateTime;

 begin
  TTT00 := Time;

  if (ImprimirIteraciones or ImprimirNumIteraciones)
                and not ArchivoSalidaAbierto then begin
    assign(F,ArchivoDeSalida);
    rewrite(F);
    ArchivoSalidaAbierto := true;
    end;

  if ImprimirIteraciones and PuedeImprimir then begin
    Writeln(F);
    Writeln(F,'Cálculo estático de flujo y K efectivo');
    Writeln(F,'BETA = ',BETA0:5:2);
    Writeln(F);
    Writeln(F,' IT  Conv Flujo    Conv Keff      Kefectivo     Sigma       Omega');
    end;

  if TipoCalculo <> adiabatico then
    MaxIter := 630
  else
    MaxIter := MaxIteraciones;

  FI := FIA;

  EPS0 := 1.0;
//  if TipoCalculo = Adiabatico then KE := 1
//                              else
  KE := KE_ANT;
  EPS := 1.0;
  iteracion := 1;

  ACC0 := 0.0;
  for i:=1 to NY do
  for j:=1 to NX do
  for iz:=1 to NZ do
  if D[iz,i,j,1] > 0.0 then
    for grupo:=1 to ng do
      ACC0 := ACC0 + FI[iz,i,j,grupo]

  else
    for grupo:=1 to ng do begin
      FI[iz,i,j,grupo] := 0.0;
      FIA[iz,i,j,grupo] := 0.0;
      end;

  OMEGA0 := 1.0;
  OMEGA := 1.0;
  NNN := 0;

  while (iteracion <= MaxIter) and (EPS >= Precision) do begin
    ACC := 0.0;

    for iz:=1 to NZ do
    for i:=1 to NY do
    for j:=1 to NX do begin

      for grupo:=1 to ng do
        FLU[grupo] := FI[iz,i,j,grupo];

      for grupo:=1 to ng do begin
        TermIzquierdo := 0.0;

        if D[iz,i,j,1] > 0.0 then
          Diagonal := diag[iz,i,j,grupo]

        else
          Diagonal := 0.0;

        if D[iz,i,j,1] > 0.0 then begin

          if (j > 1) and (D[iz,i,j-1,1] > 0.0) then
            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,i,j,1]*FI[iz,i,j-1,grupo];

          if (j < NX) and (D[iz,i,j+1,1] > 0.0) then
            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,i,j,2]*FI[iz,i,j+1,grupo];

          if (i + j) mod 2 = 0 then
          if (i > 1) and (D[iz,i-1,j,1] > 0.0) then
            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,i,j,3]*FI[iz,i-1,j,grupo];

          if (i + j) mod 2 = 1 then
          if (i < NY) and (D[iz,i+1,j,1] > 0.0) then
            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,i,j,4]*FI[iz,i+1,j,grupo];

          if (iz > 1) and (D[iz-1,i,j,1] > 0.0) then
            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,i,j,5]*FI[iz-1,i,j,grupo];

          if (iz < NZ) and (D[iz+1,i,j,1] > 0.0) then
            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,i,j,6]*FI[iz+1,i,j,grupo];

          for grupo1:=1 to ng do
            TermIzquierdo := TermIzquierdo + espectro_estatico[grupo]*NUSIGF[iz,i,j,grupo1]*FLU[grupo1]/KE;

          for grupo1:=1 to ng do
            if grupo1 <> grupo then
            TermIzquierdo := TermIzquierdo + SIGSCATT[iz,i,j,grupo1,grupo]*FLU[grupo1];

          TermIzquierdo := TermIzquierdo/Diagonal;

          FLU[grupo] := FI[iz,i,j,grupo] + BETA0*(TermIzquierdo - FI[iz,i,j,grupo]);
          end;

        end;

      for grupo:=1 to ng do
        ACC := ACC + FLU[grupo];

      for grupo:=1 to ng do
        FI[iz,i,j,grupo] := FLU[grupo];

      end;

    KE := KE*ACC/ACC0;
    EPSKE := abs(KE-KE_ANT)/KE;
    KE_ANT := KE;
    ACC0 := ACC;
    EPS := 0.0;

    NP := 0;

    for i:=1 to NY do
    for j:=1 to NX do
    for iz:=NZMIN to NZMAX do
    if D[iz,i,j,1] > 0.0 then
      for grupo:=1 to ng do begin
      ex := abs((FI[iz,i,j,grupo] - FIA[iz,i,j,grupo])/FI[iz,i,j,grupo]);
      EPS := EPS + ex;
      inc(NP);
      end;

    EPS := EPS/NP;

    if OMEGA > 0.0 then
      EPSOMEGA := abs(OMEGA - OMEGA0)/OMEGA
    else
      EPSOMEGA := 1.0;

    OMEGA0 := OMEGA;

    if (OMEGA > 1.0) and (OMEGA0 > 1.0)
        and (EPSOMEGA < 0.01) then inc(NNN)
    else NNN := 0;

    if NNN >= 3 then begin
      for i:=1 to NY do
      for j:=1 to NX do
      for iz:=NZMIN to NZMAX do
      for grupo:=1 to ng do begin
        FI[iz,i,j,grupo] := FIA[iz,i,j,grupo] + OMEGA*(FI[iz,i,j,grupo] - FIA[iz,i,j,grupo]);
        end;

      OMEGAX := OMEGA;
      NNN := 0;
      end

    else OMEGAX := 1.0;

    ALFA := EPS/EPS0;
    if ALFA < 0.999999 then OMEGA := 1.0/(1.0 - ALFA)
                       else OMEGA := 1.0;

    if PuedeImprimir then
    if ImprimirIteraciones  then
      writeln  (F,iteracion:3,'  ',
         FloatToStrF(EPS,ffexponent,5,2),'    ',
         FloatToStrF(EPSKE,ffexponent,5,2),'    ',
             KE:10:6,'  ',
            ALFA:9:5,'  ',OMEGAX:9:3);

    EPS0 := EPS;
    FIA := FI;

    inc(iteracion);
    end;

  NumIteraciones := iteracion - 1;
  UltimaConvergencia := EPS;
  PrecisionKeff := EPSKE;
  KE_ANT := KE;
  RO := (KE-1)/KE;
  if TipoCalculo = estatico then begin
    KE0 := KE;
    RO := 0.0;
    end;

  TTT00 := (Time - TTT00)*86400.0;

//  if PuedeImprimir then
  if ImprimirIteraciones or ImprimirNumIteraciones then begin
    writeln (F);

    if TipoCalculo = adiabatico then
      writeln (F,'Método Adiabático:')
    else
      writeln (F,'Cálculo Estático de Flujo:');
    writeln  (F,'Tiempo para iteraciones: ',TTT00:1:4);
    writeln (F,'Num Iteraciones: ', NumIteraciones);
    writeln (F);
    end;

  end;  (*  iteraciones  *)

(*******************************************************************************)
(*******************************************************************************)
(*******************************************************************************)
 {procedure ProdMatVect (var Res, f: OneDimension; resto, TermIndependiente: boolean);
 var
   iteracion, NNN,i,j,iz, n, NP, grupo, grupo1: integer;
   MaxIter: integer;
   Diagonal, TermIzquierdo: real;
   FLU: array[1..ng] of real;
   ACCF, ACCA, VOL, ex, ey: real;
   Result,FI: PAllLatticeMatrix;

 begin
    Result := addr(Res);
    FI := addr(f);

    VOL := VolCeldilla;
    ACCF := 0.0;
    ACCA := 0.0;

    for i:=1 to NY do
    for j:=1 to NX do
    for iz:=1 to NZ do begin

      for grupo:=1 to ng do
        FLU[grupo] := FI[iz,i,j,grupo];

      for grupo:=1 to ng do begin
        TermIzquierdo := 0.0;

        if D[iz,i,j,1] > 0.0 then
          Diagonal := diag[iz,i,j,grupo]

        else begin
          Diagonal := 0.0;
          FI[iz,i,j,grupo] := 0.0;
          Result[iz,i,j,grupo] := 0.0;
          end;

        if D[iz,i,j,1] > 0.0 then begin

          if (j > 1) and (D[iz,i,j-1,1] > 0.0) then
            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,i,j,1]*FI[iz,i,j-1,grupo];

          if (j < NX) and (D[iz,i,j+1,1] > 0.0) then
            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,i,j,2]*FI[iz,i,j+1,grupo];

          if (i + j) mod 2 = 0 then
          if (i > 1) and (D[iz,i-1,j,1] > 0.0) then
            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,i,j,3]*FI[iz,i-1,j,grupo];

          if (i < NY) and (D[iz,i+1,j,1] > 0.0) then
            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,i,j,4]*FI[iz,i+1,j,grupo];

          if (iz > 1) and (D[iz-1,i,j,1] > 0.0) then
            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,i,j,5]*FI[iz-1,i,j,grupo];

          if (iz < NZ) and (D[iz+1,i,j,1] > 0.0) then
            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,i,j,6]*FI[iz+1,i,j,grupo];

          ex := TermIzquierdo;

          for grupo1:=1 to ng do
            TermIzquierdo := TermIzquierdo + espectro_dinamico[grupo]*NUSIGF[iz,i,j,grupo1]*FLU[grupo1];

          for grupo1:=1 to ng do
            if grupo1 <> grupo then
              TermIzquierdo := TermIzquierdo + SIGSCATT[iz,i,j,grupo1,grupo]*FLU[grupo1];

          ey := (Diagonal + LS0[grupo]*VOL);

          Result[iz,i,j,grupo]:= TermIzquierdo/ey - FLU[grupo];

          if resto then
            Result[iz,i,j,grupo]:= Result[iz,i,j,grupo] + LeftSide[iz,i,j,grupo]/ey
          else if TermIndependiente then
            Result[iz,i,j,grupo]:= LeftSide[iz,i,j,grupo]/ey
          else
            Result[iz,i,j,grupo]:= -Result[iz,i,j,grupo];

      (*
          Result[iz,i,j,grupo]:= TermIzquierdo - FLU[grupo]*(Diagonal + LS0[grupo]*VOL);

          if resto then
            Result[iz,i,j,grupo]:= Result[iz,i,j,grupo] + LeftSide[iz,i,j,grupo]
          else
            Result[iz,i,j,grupo]:= -Result[iz,i,j,grupo];

     *)
          if resto then begin
//            FLU[grupo] := (TermIzquierdo + LeftSide[iz,i,j,grupo])/(Diagonal + LS0[grupo]*VOL);
            for grupo1:=1 to ng do
              ACCF := espectro_estatico[grupo]*NUSIGF0[iz,i,j,grupo1]*FLU[grupo1]*VOL + ACCF;

            ACCA := - ex + diagonal*FLU[grupo] + ACCA;

            for grupo1:=1 to ng do
              ACCA := - SIGSCATT[iz,i,j,grupo1,grupo]*FLU[grupo1]*ord(grupo1 <> grupo)  + ACCA;

            KE := ACCF/ACCA;
            end;

          end;  //  D

        end;  // grupo

//      for grupo:=1 to ng do
//        FI[iz,i,j,grupo] := FLU[grupo];

      end; //iz

  end;  }

// *****************************************************************************
 function BuscarValor (LL: TINT11; NN: integer): integer;
 var
   i: integer;

 begin
   Result := 0;
   for i:=1 to 11 do
     if LL[i] = NN then begin
       Result := i;
       end;

   end;
  {
(*******************************************************************************)
 procedure FormarILU;
 var
   i,j,iz, grupo, grupo1,LC: integer;
   Diagonal: real;
   VOL, ey: real;

// -----------------------------------------------------------------------------
 function LocCoef (iz,i,j,g: integer): integer;
 begin
   Result := (iz-1)*ng*NX*NY + (i-1)*ng*NX + (j-1)*ng + g;
   end;

// -----------------------------------------------------------------------------
{
 procedure Triangularizar;
 var
   i,j,k,kkj, kki, kk, iij, kik, kij: integer;
   ex: real;

 begin

(*   for k:=1 to NPTotal do begin
       for j:=k+1 to NPTotal do
         u[k,j] := u[k,j]/u[k,k];
       for i:=k+1 to NPTotal do
       for j:=k+1 to NPTotal do
         u[i,j] := u[i,j] - u[i,k]*u[k,j];
       end; *)

   for k:=1 to NPTotal do begin
     kk := BuscarValor(Lugar^[k],k);

     if kk > 0 then begin
       ex := valor^[k][kk];

//     for j:=k+1 to NPTotal do begin
       for kkj:=1 to 11 do begin
         j := Lugar^[k][kkj];
         if j > k then
           valor^[k][kkj] := valor^[k][kkj]/ex;
         end;

//       for j:=k+1 to NPTotal do
//       for i:=k+1 to NPTotal do
//         u[i,j] := u[i,j] - u[i,k]*u[k,j];

       for kkj:=1 to 11 do begin
         j := Lugar^[k][kkj];
         if j > k then
         for kki:=1 to 11 do begin
           i := Lugar^[k][kki];
           if i > k then begin
             kik := BuscarValor(Lugar[i],k);
             if kik > 0 then begin
               kij := BuscarValor(Lugar[i],j);
               if kij > 0 then
               valor^[i][kij] := valor^[i][kij] - valor^[i][kik]*valor^[k][kkj];
               end;
             end;
           end;
         end;

       end;
     end;

  (*
   for i:=1 to ND do
   for j:=1 to ND do
     if j < i then
       u[i,j] := 0.0
     else if A[i,j] = 0.0 then
       u[i,j] := 0.0;
  *)

     for i:=1 to NPTotal do
     for iij:=1 to 11 do begin
       j := Lugar^[i][iij];
       if j < i then
         valor^[i][iij] := 0.0;
       end;

   end;
   }
// -----------------------------------------------------------------------------
 begin
    for i:=1 to NPTotal do
    for j:=1 to 11 do begin
      Lugar^[i,j] := 0;
      valor^[i,j] := 0;
      end;

    VOL := VolCeldilla;
    LC := 0;

    for iz:=1 to NZ do
    for i:=1 to NY do
    for j:=1 to NX do begin

      for grupo:=1 to ng do begin
        inc(LC);

        if D[iz,i,j,1] > 0.0 then
          Diagonal := diag[iz,i,j,grupo]

        else begin
          Diagonal := 0.0;
          end;

        if D[iz,i,j,1] > 0.0 then begin
          ey := 1.0/(Diagonal + LS0[grupo]*VOL);

          if (j > 1) and (D[iz,i,j-1,1] > 0.0) then begin
            lugar^ [LC,6] := LocCoef(iz,i,j-1,grupo);
            valor^ [LC,6] := -coeff[grupo,iz,i,j,1]*ey;
//            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,i,j,1]*FI[iz,i,j-1,grupo];
            end;

          if (j < NX) and (D[iz,i,j+1,1] > 0.0) then begin
            lugar^ [LC,7] := LocCoef(iz,i,j+1,grupo);
            valor^ [LC,7] := -coeff[grupo,iz,i,j,2]*ey;
//            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,i,j,2]*FI[iz,i,j+1,grupo];
            end;


          if (i + j) mod 2 = 0 then
          if (i > 1) and (D[iz,i-1,j,1] > 0.0) then begin
            lugar^ [LC,8] := LocCoef(iz,i-1,j,grupo);
            valor^ [LC,8] := -coeff[grupo,iz,i,j,3]*ey;
//            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,i,j,3]*FI[iz,i-1,j,grupo];
            end;

          if (i < NY) and (D[iz,i+1,j,1] > 0.0) then begin
            lugar^ [LC,9] := LocCoef(iz,i+1,j,grupo);
            valor^ [LC,9] := -coeff[grupo,iz,i,j,4]*ey;
//            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,i,j,4]*FI[iz,i+1,j,grupo];
             end;

          if (iz > 1) and (D[iz-1,i,j,1] > 0.0) then begin
            lugar^ [LC,10] := LocCoef(iz-1,i,j,grupo);
            valor^ [LC,10] := -coeff[grupo,iz,i,j,5]*ey;
//            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,i,j,5]*FI[iz-1,i,j,grupo];
            end;

          if (iz < NZ) and (D[iz+1,i,j,1] > 0.0) then begin
            lugar^ [LC,11] := LocCoef(iz+1,i,j,grupo);
            valor^ [LC,11] := -coeff[grupo,iz,i,j,6]*ey;
//            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,i,j,6]*FI[iz+1,i,j,grupo];
            end;

          for grupo1:=1 to ng do begin
            lugar^ [LC,grupo1] := LocCoef(iz,i,j,grupo1);
            valor^ [LC,grupo1] := -espectro_dinamico[grupo]*NUSIGF[iz,i,j,grupo1]*ey;
            if grupo1 <> grupo then
              valor^ [LC,grupo1] := valor^ [LC,grupo1] - SIGSCATT[iz,i,j,grupo1,grupo]*ey;
//            TermIzquierdo := TermIzquierdo + espectro_dinamico[grupo]*NUSIGF[iz,i,j,grupo1]*FLU[grupo1];
            end;

//          for grupo1:=1 to ng do
//            if grupo1 <> grupo then
//              TermIzquierdo := TermIzquierdo + SIGSCATT[iz,i,j,grupo1,grupo]*FLU[grupo1];

//          ey := (Diagonal + LS0[grupo]*VOL);

            lugar^ [LC,grupo] := LocCoef(iz,i,j,grupo);
            valor^ [LC,grupo] := valor^ [LC,grupo] + 1.0;

//          Result[iz,i,j,grupo]:= TermIzquierdo/ey - FLU[grupo];

//          if resto then
//            Result[iz,i,j,grupo]:= Result[iz,i,j,grupo] + LeftSide[iz,i,j,grupo]/ey
//          else
//          Result[iz,i,j,grupo]:= -Result[iz,i,j,grupo];

          end;  //  D

        end;  // grupo

      end; //iz

    triangularizar;

  end;
      }
// *****************************************************************************
{ function BackRecursion (var p: OneDimension): OneDimension;
 var
   i,j,ki,kj: integer;
   ex: real;

 begin
(*      y[ND] := p[ND]/AP[ND,ND];
      for i:=ND-1 downto 1 do begin
        y[i] := p[i];
        for j:=i+1 ot NPTotal do
           y[i] := (y[i] - AP[i,j]*p[j]);
        y[i] := y[i]/AP[i,i];
        end;
  *)
//   ki := BuscarValor (lugar[NPTotal],NPTotal);
//   Result[NPTotal] := p[NPTotal]/valor[NPTotal][ki];

   for i:=NPTotal downto 1 do begin
     ki := BuscarValor ( (lugar[i]),i);

     if ki > 0 then begin
       ex := valor[i][ki];
       Result[i] := p[i];

       for kj:=1 to 11 do begin
         j := lugar[i][kj];
         if j > i then
            Result[i] := Result[i] - valor[i][kj]*p[i];
         end;

       Result[i] := Result[i]/ex;
       end

     else
       Result[i] := 0.0;

     end;

   end;
 }
(*******************************************************************************)
 procedure ProdEscalar (var R: real; var R1,R2: OneDimension);
 var
   i: integer;

 begin
   R := 0.0;
   for i:=1 to NPTotal do
     R := R + R1[i]*R2[i];

  end;

(*******************************************************************************)
 function modulo (var R: OneDimension): real;
 var
   i: integer;

 begin
   Result := 0.0;
   for i:=1 to NPTotal do
     Result := Result + R[i]*R[i];

   Result := sqrt(result);
  end;

//**************************************************************************/
//*  Resuelve el problema de cálculo del flujo para la matriz hallada      */
//*  en estado no estacionario para problemas cinético espaciales          */
//**************************************************************************/

 procedure CalculoCineticoEspacial(precision: real);
 var
   iteracion, NNN, n, i,j,iz, NP, grupo, grupo1: integer;
   MaxIter: integer;
   EPS0: real;
   Diagonal, TermIzquierdo: real;
   FLU: array[1..ng] of real;
   OMEGA, EPSOMEGA, OMEGA0, OMEGAX: real;
   KE, EPS, ALFA, ACCF, ACCA, VOL, ex: real;
   TTT00: TDateTime;
   xa,r0: POneDimension;

 begin
  TTT00 := Time;
  VOL := VolCeldilla;
//  ImprimirIteraciones := true;

  if (ImprimirIteraciones or ImprimirNumIteraciones)
                 and not ArchivoSalidaAbierto then begin
    assign(F,ArchivoDeSalida);
    rewrite(F);
    ArchivoSalidaAbierto := true;
    end;

  if ImprimirIteraciones and PuedeImprimir then begin
    Writeln(F);
    Writeln(F,'Cálculo cinético espacial por sobrerrelajación');
    Writeln(F,'BETA = ',BETA0:5:2);
    Writeln(F);
    Writeln(F,' IT  Conv Flujo  Kefectivo  Sigma    Omega');
    end;
(*
  for iz:=1 to NZ do
  for i:=1 to NY do
  for j:=1 to NX do
  for grupo:=1 to ng do begin
    FI[iz,i,j,grupo] := 1.0;
    FIA[iz,i,j,grupo] := 1.0;
    end;
  *)
  FI := FIA;
  EPS0 := 1.0;
  KE := 1.0;
  EPS := 1.0;
  iteracion := 1;
  OMEGA0 := 1.0;
  OMEGA := 1.0;
  NNN := 0;

  if DeltaT > 0.3 then MaxIter := 3000
                  else MaxIter := MaxIteraciones;

  while (iteracion <= MaxIter) and (EPS >= Precision) do begin
    ACCF := 0.0;
    ACCA := 0.0;

    for i:=1 to NY do
    for j:=1 to NX do
    for iz:=1 to NZ do begin

      for grupo:=1 to ng do
        FLU[grupo] := FI[iz,i,j,grupo];

      for grupo:=1 to ng do begin
        TermIzquierdo := 0.0;

        if D[iz,i,j,1] > 0.0 then
          Diagonal := diag[iz,i,j,grupo]

        else begin
          Diagonal := 0.0;
          FI[iz,i,j,grupo] := 0.0;
          FIA[iz,i,j,grupo] := 0.0;
          end;

        if D[iz,i,j,1] > 0.0 then begin

          if (j > 1) and (D[iz,i,j-1,1] > 0.0) then
            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,i,j,1]*FI[iz,i,j-1,grupo];

          if (j < NX) and (D[iz,i,j+1,1] > 0.0) then
            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,i,j,2]*FI[iz,i,j+1,grupo];

          if (i + j) mod 2 = 0 then
          if (i > 1) and (D[iz,i-1,j,1] > 0.0) then
            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,i,j,3]*FI[iz,i-1,j,grupo];

          if (i < NY) and (D[iz,i+1,j,1] > 0.0) then
            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,i,j,4]*FI[iz,i+1,j,grupo];

          if (iz > 1) and (D[iz-1,i,j,1] > 0.0) then
            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,i,j,5]*FI[iz-1,i,j,grupo];

          if (iz < NZ) and (D[iz+1,i,j,1] > 0.0) then
            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,i,j,6]*FI[iz+1,i,j,grupo];

          ex := TermIzquierdo;

          for grupo1:=1 to ng do
            TermIzquierdo := TermIzquierdo + espectro_dinamico[grupo]*NUSIGF[iz,i,j,grupo1]*FLU[grupo1];

          for grupo1:=1 to ng do
            if grupo1 <> grupo then
              TermIzquierdo := TermIzquierdo + SIGSCATT[iz,i,j,grupo1,grupo]*FLU[grupo1];

          FLU[grupo] := (TermIzquierdo + LeftSide[iz,i,j,grupo])/(Diagonal + LS0[grupo]*VOL);

          FLU[grupo] := FIA[iz,i,j,grupo] + BETA0*(FLU[grupo] - FIA[iz,i,j,grupo]);

          for grupo1:=1 to ng do
            ACCF := espectro_estatico[grupo]*NUSIGF0[iz,i,j,grupo1]*FLU[grupo1]*VOL + ACCF;

          ACCA := - ex + diagonal*FLU[grupo] + ACCA;

          for grupo1:=1 to ng do
            ACCA := - SIGSCATT[iz,i,j,grupo1,grupo]*FLU[grupo1]*ord(grupo1 <> grupo)  + ACCA;
          end;  //  D

        end;  // grupo

      for grupo:=1 to ng do
        FI[iz,i,j,grupo] := FLU[grupo];

      end; //iz

    KE := ACCF/ACCA;
    EPS := 0.0;
    NP := 0;

    for iz:=1 to NZ do
    for i:=1 to NY do
    for j:=1 to NX do
    if D[iz,i,j,1] > 0.0 then begin
      inc(NP);
      for grupo:=1 to ng do
        EPS := EPS + abs((FI[iz,i,j,grupo] - FIA[iz,i,j,grupo])/FI[iz,i,j,grupo])
      end;

    EPS := EPS/NP;

    if OMEGA > 0.0 then
      EPSOMEGA := abs(OMEGA - OMEGA0)/OMEGA
    else
      EPSOMEGA := 1.0;

    OMEGA0 := OMEGA;

    if (OMEGA > 1.0) and (OMEGA0 > 1.0)
        and (EPSOMEGA < 0.01) then inc(NNN)
    else NNN := 0;

    if NNN >= 3 then begin
      for i:=1 to NY do
      for j:=1 to NX do
      for iz:=1 to NZ do
      for grupo:=1 to ng do begin
        FI[iz,i,j,grupo] := FIA[iz,i,j,grupo] + OMEGA*(FI[iz,i,j,grupo] - FIA[iz,i,j,grupo]);
        end;

      OMEGAX := OMEGA;
      NNN := 0;
      end

    else OMEGAX := 1.0;

    ALFA := EPS/EPS0;
    if ALFA < 0.999999 then OMEGA := 1.0/(1.0 - ALFA)
                       else OMEGA := 1.0;

    if PuedeImprimir then
    if ImprimirIteraciones then
      writeln  (F,iteracion:3,'  ',FloatToStrF(EPS,ffexponent,5,2),'  ',KE:1:6,
                     '  ',ALFA:1:5,'  ',OMEGAX:6:3);

    EPS0 := EPS;
    FIA := FI;

    inc(iteracion);
    end;  // while

  TTT00 := (Time - TTT00)*86400.0;

  NumIteraciones := iteracion - 1;
  UltimaConvergencia := EPS;

  if PuedeImprimir then
  if ImprimirIteraciones or ImprimirNumIteraciones then begin
    new (xa);
    new (r0);
    writeln (F);
    move(FI,xa^,SizeOf(FI)); // xa^ := FI;
    ProdMatVect (r0^,xa^,true,false);
    writeln (F,'Método Directo:');
    writeln (F,'Error del Resto: ',FloatToStrF(modulo(r0^)/modulo(xa^),ffexponent,5,2));

    writeln  (F,'Tiempo para iteraciones: ',TTT00:1:3);
    writeln (F,'Num Iteraciones: ', NumIteraciones);
    writeln (F);
    dispose (xa);
    dispose (r0);
    end;

  RO := (KE-1)/KE;

  end; (*  CalculoCineticoEspacial  *)

(*******************************************************************************)
 procedure precursores(DT: real);
 var
   iz,i,j,prec,grupo: integer;
   FF: real;

 begin
   for iz:=1 to NZ do
   for i:=1 to NY do
   for j :=1 to NX do
   if D[iz,i,j,1] > 0.0 then begin
     FF := 0.0;
     for grupo:=1 to ng do
       FF := FF + NUSIGF0[iz,i,j,grupo]*FIANT[iz,i,j,grupo];

     for prec:=1 to NPREC do
     if DT = 0.0 then
       ConcPrecursor[prec,iz,i,j] := FF*BetaLambda[prec]/KE0
     else
       ConcPrecursor[prec,iz,i,j] :=
         (ConcPrecursor[prec,iz,i,j] + FF*BetaNuclear[prec]*DT)/(1 + Lambda[prec]*DT);
     end;

  end;  (*  precursores  *)

(**********************************************************)
procedure linear_variation(var cn: real; var con: PrecursorType;
         var der_log, fuen: real; r0,r1,dt: real);
var
  i,np,nn: integer;
  du,dro,ro,cna,bet,der_logd,fact,eps,tt: real;
  cona: PrecursorType;
  coef: PrecursorType;

(* --------------------------------------------------------- *)
(*  Calculation of one step supposing that reactivity varies *)
(*     linearly from ro up to r1 during interval du          *)
(* --------------------------------------------------------- *)

 procedure OneStep (var ro,r1,du: real);
 label fuera;
 var
   cnd,cnd1,cnd2,term,finv,ddt: real;
   i,k,l: integer;

 begin
   der_logd := 0;

   cnd := cna;
   cnd1 := cna;
   term := cna;
   finv := 1;
   ddt := 1;
   k := 0;

 (* Using a Taylor development up to order 30. It is expected  *)
 (* to converge in no more than 6 or 7 iterations because of   *)
 (* the criteria to choose the small step du                   *)

   while k <= 30 do begin
     if not (abs(term/cna) > eps) then goto fuera;
     inc(k);
     cnd2 := cnd1;
     cnd1 := cnd;
     cnd := (ro-bet)/MeanLifeTime*cnd1;

     if k = 1 then cnd := cnd + fuen;

     for l:=1 to nprec do
       cnd := cnd + lambda[l]*coef[l];

     cnd := cnd + (k-1)*r1/MeanLifeTime*cnd2;

     for i:=1 to NPREC do
       coef[i] := -lambda[i]*coef[i] + BetaNuclear[i]/MeanLifeTime*cnd1;

     der_logd := der_logd + cnd*finv*ddt;
     finv := finv/k;
     ddt := ddt*du;
     term := cnd*finv*ddt;
     cna := cna + term;

     for i:=1 to NPREC do
       cona[i] := cona[i] + coef[i]*finv*ddt;

     end;

 fuera:
//   if k > 30 then error
//    ('Point kinetics with linear vatiation of reactivity does not work');

   end; (* OneStep *)

(* ------------------------------------- *)
begin        (*   linear_variation  *)
  fact := 0.05;
  eps := 1.0e-8;

  bet := BetaNuclearTotal;

 (* Empyrical formula for calculating a small step
    for resolution of the point kinetic equations
    for linear reactivity variation *)
  du := MeanLifeTime/((abs(r0)+abs(r1))+bet)*fact;
  np := trunc(dt/du)+1;
  du := dt/np;
  dro := r1*dt/np;
  ro := r0-dro;

  for i:=1 to NPREC do cona[i] := con[i];
  cna := cn;
  tt := 0;

  for nn:=1 to np do begin
    tt := tt + du;
    ro := ro + dro;
    for i:=1 to NPREC do coef[i] := cona[i];

  (* Calculation of one step supposing that reactivity varies
     linearly from ro up to r1 during interval du *)
    OneStep (ro,r1,du);
    end;

  for i:=1 to NPREC do con[i] := cona[i];
  cn := cna;
  der_log := der_logd/cna;

//  cona := nil;
//  coef := nil;
  end;  (* linear_variation *)

(**********************************************************)
procedure Backwards_differences
    (var cn: real; var con: PrecursorType; var der_log: real; RO,fuen,dt: real);
var
  acc0,acc: real;
  i: integer;

begin
  acc0 := 0.0;
  acc := 0.0;
  der_log := cn;

  for i:=1 to 6 do begin
    acc0 := acc0 + lambda[i]*con[i]/(1 + Lambda[i]*DT);
    acc := acc + BetaNuclear[i]/(1 + Lambda[i]*DT);
    end;

  cn := (cn + dt*(acc0 + fuen))/(1.0 - dt/MeanLifeTime*(RO - acc));

  for i:=1 to 6 do
    con[i] := (con[i] + dt*BetaNuclear[i]*cn/MeanLifeTime)/(1.0 + lambda[i]*dt);

  der_log := (cn - der_log)/((cn+der_log)*0.5*dt);
  end;

(*******************************************************************************)
 procedure AdiabaticMethod(DT,DTXenon: real);
 var
   prec,iz,i,j, grupo: integer;
   ACC, ACCF, der_log, cn, fuen, r0, r1: real;
   conc: PrecursorType;
//   concx: PrecursorType;

 begin
   InterpolateXS;

   r0 := RO;
   fuen := 0.0;
   ACC := 0.0;
   ACCF := 0.0;
   for prec:=1 to NPREC do conc[prec] := 0.0;

   for iz:=NZMIN to NZMAX do
   for i:=1 to NY do
   for j:=1 to NX do
   if D[iz,i,j,1] > 0.0 then begin
     for grupo:=1 to ng do begin
       ACC := ACC + FI[iz,i,j,grupo]/VEL[grupo];
       ACCF := ACCF + FI[iz,i,j,grupo]*NUSIGF0[iz,i,j,grupo]/KE0;
       end;

     for prec:=1 to NPREC do
       conc[prec] := conc[prec] + ConcPrecursor[prec,iz,i,j];

     end;

   MeanLifeTime := ACC/ACCF;
   cn := ACC;
{
   for prec:=1 to NPREC do
      concx[prec] := BetaNuclear[prec]*ACC/(lambda[prec]*MeanLifeTime);

//   conc := concx;

  if ImprimirIteraciones and not ArchivoSalidaAbierto then begin
    assign(F,ArchivoDeSalida);
    rewrite(F);
    ArchivoSalidaAbierto := true;
    end;

   writeln(F);
   writeln(F);
//   writeln (F,'Const Puntuales T = ',TT:1:4);
   writeln (F,' L = ',MeanLifeTime);
   writeln (F,'cn = ',cn);
   writeln (F,conc[1],' ',conc[2],' ',conc[3],' ',conc[4],' ',conc[5],' ',conc[6]);
   writeln(F);
   writeln(F);
   writeln (F,concx[1],' ',concx[2],' ',concx[3],' ',concx[4],' ',concx[5],' ',concx[6]);

//   conc := concx;
   CloseFile(F);
   halt(0);                       }
   FIANT := FI;
   MeshContants (adiabatico,DT);
   iteraciones(adiabatico,PrecisionCinetico);

   r1 := RO;

//   if variacion_lineal then
//     linear_variation(cn, conc, der_log, fuen, r0,r1,DT)
//   else
   Backwards_differences (cn, conc, der_log, RO, 1.0E6, DT);

   DerivadaLogaritmica := der_log;

   ACC := 0.0;

   for iz:=1 to NZ do
   for i:=1 to NY do
   for j:=1 to NX do
   if D[iz,i,j,1] > 0.0 then
     for grupo:=1 to ng do
       ACC := ACC + FI[iz,i,j,grupo]/VEL[grupo];

   ACC := cn/ACC;

   for iz:=1 to NZ do
   for i:=1 to NY do
   for j:=1 to NX do
   if D[iz,i,j,1] > 0.0 then
     for grupo:=1 to ng do
       FI[iz,i,j,grupo] := FI[iz,i,j,grupo]*ACC;

   potencias(0.0,DT);
   precursores(DT);
   if ConXenon then XENON(DTXenon);
   if ConSamario then SAMARIO(DTXenon);
   end;

(*******************************************************************************)
 procedure DirectMethod(DT, DTXenon: real);
 var
   iz,i,j, grupo: integer;
   CONC1,CONC2: real;

 begin
   CONC1 := 0.0;

   for iz:=1 to NZ do
   for i:=1 to NY do
   for j:=1 to NX do
   if D[iz,i,j,1] > 0.0 then
     for grupo:=1 to ng do
       CONC1 := CONC1 + FI[iz,i,j,grupo]/VEL[grupo];

   FIANT := FI;
   InterpolateXS;
   MeshContants (directo,DT);
   CalculoCineticoEspacial(PrecisionCinetico);
   potencias(0.0,DT);
   precursores(DT);

   CONC2 := 0.0;

   for iz:=1 to NZ do
   for i:=1 to NY do
   for j:=1 to NX do
   if D[iz,i,j,1] > 0.0 then
     for grupo:=1 to ng do
       CONC2 := CONC2 + FI[iz,i,j,grupo]/VEL[grupo];

   DerivadaLogaritmica := (CONC2 - CONC1)/(CONC2*DT);

   if ConXenon then XENON(DTXenon);
   if ConSamario then SAMARIO(DTXenon);
   end;

//**************************************************************************/
//*  Evoluciona un paso DT en el tiempo para un ciclo de Xenón             */
//**************************************************************************/

 procedure XenonStaticCalculation(POWER, DTXenon: real);
 begin
   UNIDAD := 3600.0;
   InterpolateXS;
   FIANT := FI;
   MeshContants (estatico,0.0);
   iteraciones(estatico, PrecisionEstacionario);
   potencias(POWER,0.0);
   precursores(0.0);
   if ConXenon and (DTXenon > 0.0) then XENON(DTXenon);
   if ConSamario and (DTXenon > 0.0) then SAMARIO(DTXenon);
   end;

//**************************************************************************/
//*  singleiza un cálculo estático inicial                                   */
//**************************************************************************/

 procedure StaticCalculation(POWER: real);
 begin
   PuedeImprimir := true;
   InterpolateXS;
// imprimir_XS_red(1,19,1,18,7,7);

   MeshContants (estatico,0.0);
// CompararConPuma (1);
   iteraciones(estatico, PrecisionEstacionario);
   potencias(POWER,0.0);
   FIANT := FI;
   precursores(0.0);

   if ConXenon then XENON(0.0);
   if ConSamario then SAMARIO(0.0);
   PotInit := PotTotal;
   KeffInit := KE0;
   end;

(*******************************************************************************)
 procedure TransferirEstado(var TTT: real; grabar: boolean; var mem: MemType);
 var
   N, NI: integer;
   NCM, NCT, NCFM, NCICM: integer;
   NPAux, NCPrec, NAll, NTotal: integer;

 (* ------------------------------------------------------ *)
 procedure G(var V; N: integer);
 var
   PMem: ^MemType;
   i,j: integer;

 begin
   if grabar then begin
     PMem := addr(V);

     j:= 1;
     for i:=NTotal to NTotal+N-1 do begin
       Mem[i] := PMem^[j];
       inc(j);
       end;

//     BlockWrite(SF,V,N)
     end

   else begin
//     BlockRead(SF,V,N);
     PMem := addr(V);

     j:= 1;
     for i:=NTotal to NTotal+N-1 do begin
       PMem^[j] := Mem[i];
       inc(j);
       end;

     end;

   NTotal := NTotal + N;
   end;

 (* ------------------------------------------------------ *)
 begin
   NTotal := 0;

   N := SizeOf(real);
   NI := SizeOf(integer);
   NCM := SizeOf(CoreMatrix);
   NCT := SizeOf(ChannelType);
   NCFM := SizeOf(CoreFluxMatrix);
   NCICM := SizeOf(IntCoreMatrix);

   G(TTime,N);
   TTT := TTime;
   G(DeltaT,N);
   G(PotTotal,N); G(PotResidual,N);  G(PotInit,N);   G(RO,N);
   G(KeffInit,N);
   G(DerivadaLogaritmica,N); G(FForma,N);
   G(PMax,N);
   G(CanalPMax,NI); G(TrozoPMax,NI);
   G(PotMax,N);
   G(CanalPotMax,NI);

   G(ConcXEMedia,N); G(ConcXenon0,N);
   G(ConcSmMedia,N); G(ConcSamario0,N);

   G(TComb,NCM);
   G(DensRefr,NCM);
   G(TempRefr,NCM);
   G(PPMBoro,NCM);
   G(POT,NCM);
   G(potcan,NCT);
   G(flujo,NCFM);

   G(KE0,N);

   G(Barra_E5,N);
   G(Barra_E4,N);   G(Barra_D6,N);   G(Barra_F5,N);
   G(Barra_D3,N);   G(Barra_C8,N);   G(Barra_H4,N);
   G(Barra_H2,N);   G(Barra_B5,N);   G(Barra_E8,N);
   G(Barra_F2,N);   G(Barra_B7,N);   G(Barra_G6,N);
   G(Barra_F4,N);   G(Barra_D5,N);   G(Barra_E6,N);
   G(Barra_E2,N);   G(Barra_B8,N);   G(Barra_H5,N);
   G(Barra_H3,N);   G(Barra_C4,N);   G(Barra_D8,N);
   G(Barra_G2,N);   G(Barra_B6,N);   G(Barra_F7,N);
   G(Barra_D4,N);   G(Barra_D7,N);   G(Barra_G4,N);

   G(KE0,N);
   G(BurnUp,NCM);

   G(ConcXE,NCM);
   G(ConcI,NCM);

   G(ConcPt,NCM);
   G(ConcSm,NCM);

   G(NumMaterial,NCICM);

   G(NumMaterialBarra,NCICM);
   G(FraccionBarra,NCM);

   G(PotRes,NCM);
   G(PotInst,NCM);

   G(PotInstAnt,NCM);
   NPAux := SizeOf(PAux); G(PAux,NPAux);
   NCPrec := SizeOf(ConcPrecursor); G(ConcPrecursor,NCPrec);
   NAll := SizeOf(FI); G(FI,NAll); G(FIA,NAll);
   G(KE_ANT,N);
   end;

 (********************************************************************************)
  end.
