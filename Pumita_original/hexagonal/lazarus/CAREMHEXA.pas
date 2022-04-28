 (*$R+*)

unit CAREMHEXA;
interface
uses Dialogs,
  math,
  interpol,
  sysutils;

const
  NTROZOS = 14;
  NCANALES = 61;
  ng = 5;

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
 ArchivoIntMultiple: string;        // Archivo de interpolación múltiple.
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
 ImprimirPorCanal:         // Muestra la distribución de potencia por canal
       boolean;
 UNIDAD: real;

// Valores centrales de tabla
 ConcXE0: real;
 TComb0: real;
 DensRefr0: real;
 TempRefr0: real;

 PrecisionEstacionario, PrecisionCinetico: real;
 PuedeImprimir: boolean;
 FactPuma: real;

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
//  NX = 39;
//  NY = 24;
  NRegionesReactor = 61;
  NTotalRegiones = 91;

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
  LongTableLine0 = 321; 
  NumMaxLineasTabla = 49;


Type
  TableLine = array[1..LongTableLine0] of real;
  TableType = array[1..NumMaxLineasTabla] of TableLine;
  LatticeAxialLine = array[1..NZ] of real;
  AllLatticeMatrix = array[1..NZ,1..NTotalRegiones] of GroupType;
  LatticeMatrix = array[1..NZ,1..NTotalRegiones] of real;
  ScattMatrix = array[1..NZ,1..NTotalRegiones] of ScattType;
  PrecursorType = array[1..NPREC] of real;
  BetaLamdaType = array[1..NPREC,1..ng] of real;

const
  DSIGABoro = 5.0E-2*0.00055;

const
  ALTURA = 8.0;
  DELTAZ = 10.0;
  APOTEMA = 16.0/2.0;
  LADO = 9.237604307034012;
  SUP_HEXA = 36.95041722813605*6;
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


  vecino: array [1..91,1..6] of integer =
 (*  1  *)                               ((64,62,63,2,5,3),
 (*  2  *)                         (1,63,65,4,8,5),(66,64,1,5,9,6),
 (*  4  *)                  (2,65,67,7,12,8),(3,1,2,8,13,9),(68,66,3,9,14,10),
 (*  7  *)            (4,67,69,11,16,12),(5,2,4,12,17,13),(6,3,5,13,18,14),(70,68,6,14,19,15),
 (* 11  *)   (7,69,71,73,20,16),(8,4,7,16,21,17),(9,5,8,17,22,18),(10,6,9,18,23,19), (72,70,10,19,24,74),
 (* 16  *)       (12,7,11,20,25,21),(13,8,12,21,26,22),(14,9,13,22,27,23),(15,10,14,23,28,24),
 (* 20  *)   (16,11,73,75,29,25),(17,12,16,25,30,26),(18,13,17,26,31,27),(19,14,18,27,32,28),(74,15,19,28,33,76),
 (* 25  *)       (21,16,20,29,34,30),(22,17,21,30,35,31),(23,18,22,31,36,32),(24,19,23,32,37,33),
 (* 29  *)   (25,20,75,77,38,34),(26,21,25,34,39,35),(27,22,26,35,40,36),(28,23,27,36,41,37),(76,24,28,37,42,78),
 (* 34  *)       (30,25,29,38,43,39),(31,26,30,39,44,40),(32,27,31,40,45,41),(33,28,32,41,46,42),
 (* 38  *)   (34,29,77,79,47,43),(35,30,34,43,48,44),(36,31,35,44,49,45),(37,32,36,45,50,46),(78,33,37,46,51,80),
 (* 43  *)       (39,34,38,47,52,48),(40,35,39,48,53,49),(41,36,40,49,54,50),(42,37,41,50,55,51),
 (* 47  *)   (43,38,79,81,83,52),(44,39,43,52,56,53),(45,40,44,53,57,54),(46,41,45,54,58,55),(80,42,46,55,84,82),
 (* 52  *)       (48,43,47,83,85,56),(49,44,48,56,59,57),(50,45,49,57,60,58),(51,46,50,58,86,84),
 (* 56  *)              (53,48,52,85,87,59),(54,49,53,59,61,60),(55,50,54,60,88,86),
 (* 59  *)                      (57,53,56,87,89,61),(58,54,57,61,90,88),
 (* 61  *)                                 (60,57,59,89,91,90),
 (* 62  *)             (0,0,0,63,1,64),
 (* 63  *)   (62,0,0,65,2,1),(0,0,62,1,3,66),
 (* 65  *)   (63,0,0,67,4,2),(0,0,64,3,6,68),
 (* 67  *)   (65,0,0,69,7,4),(0,0,66,6,10,70),
 (* 69  *)   (67,0,0,71,11,7),(0,0,68,10,15,72),
 (* 71  *)   (69,0,0,0,73,11),(0,0,70,15,74,0),
 (* 73  *)   (11,71,0,0,75,20),(0,72,15,24,76,0),
 (* 75  *)   (20,73,0,0,77,29),(0,74,24,33,78,0),
 (* 77  *)   (29,75,0,0,79,38),(0,76,33,42,80,0),
 (* 79  *)   (38,77,0,0,81,47),(0,78,42,51,82,0),
 (* 81  *)   (47,79,0,0,0,83), (0,80,51,84,0,0),
 (* 83  *)   (52,47,81,0,0,85),(82,51,55,86,0,0),
 (* 85  *)   (56,52,83,0,0,87),(84,55,58,88,0,0),
 (* 87  *)   (59,56,85,0,0,89),(86,58,60,90,0,0),
 (* 89  *)   (61,59,87,0,0,91),(88,60,61,91,0,0),
 (* 91  *)            (90,61,89,0,0,0));


  reticula: array [1..141,1..4] of integer =

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

// *** FILA 3
// ***  region 62
   (18, 23,  3, 4),

// *** FILA 4
// ***  region 63
   (15, 19,  4, 5),
// ***  region 64
   (22, 26,  4, 5),

// *** FILA 5
// ***  region 65
   (12, 16,  5, 6),
// ***  region 66
   (25, 29,  5, 6),

// *** FILA 6
// ***  region 67
   (9, 13,  6, 7),
// ***  region 68
   (28, 32,  6, 7),

// *** FILA 7
// ***  region 69
   (6, 10,  7, 8),
// ***  region 70
   (31,35,  7, 8),

// *** FILA 8
// ***  region 71
   (5,  7,  8, 9),
// ***  region 72
   (34,36,  8, 9),

// *** FILA 9
// ***  region 73
   (5,  7,  9,10),
// ***  region 74
   (34,36,  9,10),

// *** FILA 10
// ***  region 75
   (5,  7, 10,11),
// ***  region 76
   (34,36, 10,11),

// *** FILA 11
// ***  region 77
   (5,  7, 11,12),
// ***  region 78
   (34,36, 11,12),

// *** FILA 12
// ***  region 79
   (5,  7, 12,13),
// ***  region 80
   (34,36, 12,13),

// *** FILA 13
// ***  region 81
   (5,  7, 13,14),
// ***  region 82
   (34,36, 13,14),

// *** FILA 14
// ***  region 83
   (5,  7, 14,15),
// ***  region 84
   (34,36, 14,15),

// *** FILA 15
// ***  region 85
   (5,  7, 15,16),
// ***  region 86
   (34,36, 15,16),

// *** FILA 16
// ***  region 87
   (5,  7, 16,17),
// ***  region 88
   (34,36, 16,17),

// *** FILA 17
// ***  region 89
   (5,  7, 17,18),
// ***  region 90
   (34,36, 17,18),

// *** FILA 18
// ***  region 91
   (6, 10, 18,19),
// ***  region 92
   (31,35, 18,19),

// *** FILA 19
// ***  region 93
   (9, 13, 19,20),
// ***  region 94
   (28,32, 19,20),

// *** FILA 20
// ***  region 95
   (12,16, 20,21),
// ***  region 96
   (25,29, 20,21),

// *** FILA 21
// ***  region 97
   (15,19, 21,22),
// ***  region 98
   (22,26, 21,22),

// *** FILA 22
// ***  region 99
   (18, 23, 22,23),


// ********************************
// *****   CONTORNO EXTERNO   *****
// ********************************

// *** FILA 1 Y 2
// ***  region 100
   (1,  40, 1, 3),

// *** FILA 3
// ***  region 101
   (1,  18,  3, 4),
// ***  region 102
   (23, 40,  3, 4),

// *** FILA 4
// ***  region 103
   (1,  15,  4, 5),
// ***  region 104
   (26, 40,  4, 5),

// *** FILA 5
// ***  region 105
   (1,  12,  5, 6),
// ***  region 106
   (29, 40,  5, 6),

// *** FILA 6
// ***  region 107
   (1,  9,   6, 7),
// ***  region 108
   (32, 40,  6, 7),

// *** FILA 7
// ***  region 109
   (1,  6,   7, 8),
// ***  region 110
   (35, 40,  7, 8),

// *** FILA 8
// ***  region 111
   (1, 5,  8,9),
// ***  region 112
   (36,40, 8,9),

// *** FILA 9
// ***  region 113
   (1, 5,  9,10),
// ***  region 114
   (36,40, 9,10),

// *** FILA 10
// ***  region 115
   (1,  5,  10,11),
// ***  region 116
   (36, 40, 10,11),

// *** FILA 11
// ***  region 117
   (1,  5,  11,12),
// ***  region 118
   (36, 40, 11,12),

// *** FILA 12
// ***  region 119
   (1,  5, 12,13),
// ***  region 120
   (36, 40,12,13),

// *** FILA 13
// ***  region 121
   (1,  5, 13,14),
// ***  region 122
   (36, 40, 13,14),

// *** FILA 14
// ***  region 123
   (1,  5, 14,15),
// ***  region 124
   (36, 40, 14,15),

// *** FILA 15
// ***  region 125
   (1,  5, 15,16),
// ***  region 126
   (36, 40,15,16),

// *** FILA 16
// ***  region 127
   (1,  5, 16,17),
// ***  region 128
   (36, 40, 16,17),

// *** FILA 17
// ***  region 129
   (1,  5, 17,18),
// ***  region 130
   (36, 40, 17,18),

// *** FILA 18
// ***  region 131
   (1,  6, 18,19),
// ***  region 132
   (35, 40, 18,19),

// *** FILA 19
// ***  region 133
   (1,  9, 19,20),
// ***  region 134
   (32, 40, 19,20),

// *** FILA 20
// ***  region 135
   (1,  12, 20,21),
// ***  region 136
   (29, 40, 20,21),

// *** FILA 21
// ***  region 137
   (1,  15, 21,22),
// ***  region 138
   (26, 40, 21,22),

// *** FILA 22
// ***  region 139
   (1,  18, 22,23),
// ***  region 140
   (23, 40, 22,23),

// *** FILA 23 Y 24
// ***  region 141
   (1,  40, 23,25)) ;


  XS_REFL_Inf0: TipoSeccionesEficaces =
 (2.8013E+00, 1.4324E+00, 7.5381E-01, 3.0550E-01, 1.7568E-01,
  4.0668E-04, 8.1694E-02, 5.1308E-04, 0.0000E+00, 0.0000E+00,
  0.0000E+00, 8.0729E-06, 1.2447E-01, 1.2175E-05, 1.4402E-06,
  0.0000E+00, 0.0000E+00, 7.8181E-04, 1.0667E-01, 9.2764E-03,
  0.0000E+00, 0.0000E+00, 5.9855E-04, 7.1300E-03, 4.4451E-01,
  0.0000E+00, 0.0000E+00, 6.1834E-06, 4.9102E-01, 1.4114E-02,
  0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,
  0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00);

  XS_REFL_Sup0: TipoSeccionesEficaces =
 (3.2434E+00, 1.6599E+00, 8.7306E-01, 3.5640E-01, 2.0363E-01,
  3.4877E-04, 7.0647E-02, 4.4366E-04, 0.0000E+00, 0.0000E+00,
  0.0000E+00, 6.9592E-06, 1.0723E-01, 1.0488E-05, 1.2403E-06,
  0.0000E+00, 0.0000E+00, 6.7396E-04, 9.1956E-02, 7.9769E-03,
  0.0000E+00, 0.0000E+00, 6.3186E-04, 6.0791E-03, 3.6955E-01,
  0.0000E+00, 0.0000E+00, 8.7021E-06, 4.7060E-01, 1.2070E-02,
  0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,
  0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00);

  XS_REFL_Rad0: TipoSeccionesEficaces =
 (2.8013E+00, 1.4324E+00, 7.5381E-01, 3.0550E-01, 1.7568E-01,
  4.0668E-04, 8.1694E-02, 5.1308E-04, 0.0000E+00, 0.0000E+00,
  0.0000E+00, 8.0729E-06, 1.2447E-01, 1.2175E-05, 1.4402E-06,
  0.0000E+00, 0.0000E+00, 7.8181E-04, 1.0667E-01, 9.2764E-03,
  0.0000E+00, 0.0000E+00, 5.9855E-04, 7.1300E-03, 4.4451E-01,
  0.0000E+00, 0.0000E+00, 6.1834E-06, 4.9102E-01, 1.4114E-02,
  0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,
  0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00);


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
(*  ((-4.99311E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00),
  ( 1.83690E-01,-4.22252E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00),
  ( 4.31017E-02, 1.86218E-01,-2.69085E-01, 2.55694E-03, 9.57966E-04),
  ( 2.49466E-02, 5.12368E-02, 1.22192E-01,-2.55231E-01, 2.28923E-01),
  ( 1.07540E-02, 2.15486E-02, 4.16998E-02, 1.87461E-01,-2.70710E-01)); *)

   alfaaxialSup0: MatRespType =
 ((-4.99999E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00),
  ( 1.82907E-01,-4.22387E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00),
  ( 4.14578E-02, 1.86035E-01,-2.67089E-01, 3.19515E-03, 1.21206E-03),
  ( 2.15478E-02, 4.89087E-02, 1.22441E-01,-2.43686E-01, 2.47177E-01),
  ( 7.95719E-03, 1.76834E-02, 3.66747E-02, 1.74610E-01,-2.89138E-01));
(* ((-5.02749E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00),
  ( 1.82336E-01,-4.23366E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00),
  ( 4.11860E-02, 1.85696E-01,-2.67512E-01, 3.18907E-03, 1.20890E-03),
  ( 2.13364E-02, 4.86422E-02, 1.22202E-01,-2.43896E-01, 2.47080E-01),
  ( 7.87827E-03, 1.75838E-02, 3.65862E-02, 1.74541E-01,-2.89192E-01));*)

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

  Coeff: array[1..ng,1..NZ,1..NTotalRegiones,1..8] of real;

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
     ShowMessage ('Matriz singular 1');
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
       ShowMessage ('Matriz singular 2');
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
   s := SUP_HEXA;

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
  procedure MeshContants (TipoCalculo: TipoCalculoType; DT: real);
  const
    ALFA = 5.000; //0.46922;
//    ALFA = 0.41300;

  var
    canal, trozo, grupo, grupo1, i, j, k, caso,
              N00, iz, i1, i2, j1, j2, prec, nreg, borde, num_vecino: integer;
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
    for nreg:=1 to NTotalRegiones do
    if (D[iz,nreg,1] > 0.0)  and not ((iz < NZMIN) or (iz > NZMAX)) then begin
      SIGTOT[iz,nreg,grupo] :=  SIGTOT[iz,nreg,grupo]*VOL;
      NUSIGF[iz,nreg,grupo] :=  NUSIGF[iz,nreg,grupo]*VOL;
      for grupo1:=1 to ng do
        SIGSCATT[iz,nreg,grupo,grupo1] := SIGSCATT[iz,nreg,grupo,grupo1]*VOL;
      end;

    if ConMatrizDeRespuesta then
    for nreg:=1 to NTotalRegiones do
    if D[1,nreg,1] > 0.0 then begin
      hallar_contorno(D[1,nreg], alfaaxialInf0, ESCAPE_ABAJO);

      for grupo:=1 to ng do
      for grupo1:=1 to ng do
      if grupo = grupo1 then
        SIGTOT[1,nreg,grupo] := SIGTOT[1,nreg,grupo] - ESCAPE_ABAJO[grupo,grupo1]
      else
        SIGSCATT[1,nreg,grupo,grupo1] := SIGSCATT[1,nreg,grupo,grupo1]
                          +  ESCAPE_ABAJO[grupo,grupo1];

      end;

    if ConMatrizDeRespuesta then
    for nreg:=1 to NTotalRegiones do
    if D[NZ,nreg,1] > 0.0 then begin
      hallar_contorno(D[NZ,nreg], alfaaxialSup0, ESCAPE_ARRIBA);

      for grupo:=1 to ng do
      for grupo1:=1 to ng do
      if grupo = grupo1 then
        SIGTOT[NZ,nreg,grupo] := SIGTOT[NZ,nreg,grupo] - ESCAPE_ARRIBA[grupo,grupo1]
      else
        SIGSCATT[NZ,nreg,grupo,grupo1] := SIGSCATT[NZ,nreg,grupo,grupo1]
                          +  ESCAPE_ARRIBA[grupo,grupo1];
      end;

    for grupo:=1 to ng do
    for iz:=1 to NZ do
    for nreg:=1 to NTotalRegiones do
    if (D[iz,nreg,1] = 0.0)  or (iz < NZMIN) or (iz > NZMAX) then begin
      diag[iz,nreg,grupo] := 0.0;
      diag0[iz,nreg,grupo] := 0.0;

      for k:=1 to 6 do
        coeff[grupo,iz,nreg,k] := 0.0;

      end

    else begin

      diag[iz,nreg,grupo] := 0.0; //SIGTOT[grupo,iz,i,j];

      for borde:=1 to 6 do begin
        num_vecino := vecino[nreg,borde];

        if (num_vecino > 0) and (D[iz,num_vecino,1] > 0.0) then begin
          EX := SUP_LADO/(APOTEMA/D[iz,num_vecino,grupo] + APOTEMA/D[iz,nreg,grupo]);
          diag[iz,nreg,grupo] := diag[iz,nreg,grupo] + EX;
          coeff[grupo,iz,nreg,borde] := EX;
          end

        else begin
          EX := SUP_LADO*ALFA/(1.0 + APOTEMA*ALFA/D[iz,nreg,grupo]);
          coeff[grupo,iz,nreg,borde] := EX;
          diag[iz,nreg,grupo] := diag[iz,nreg,grupo] + EX;
          end;

        end;


      if (iz > NZMIN) and (D[iz-1,nreg,1] > 0.0) then begin
        EX := 2.0*SUP_HEXA/(DELTAZ/D[iz-1,nreg,grupo] + DELTAZ/D[iz,nreg,grupo]);
        diag[iz,nreg,grupo] := diag[iz,nreg,grupo] + EX;
        coeff[grupo,iz,nreg,7] := EX;
        end

      else if not ConMatrizDeRespuesta then begin
        EX := (SUP_HEXA*ALFA/(1.0 + DELTAZ*ALFA*0.5/D[iz,nreg,grupo]));
        coeff[grupo,iz,nreg,7] := EX;
        diag[iz,nreg,grupo] := diag[iz,nreg,grupo] + EX;
        end;


      if (iz < NZMAX) and (D[iz+1,nreg,1] > 0.0) then begin
        EX := 2.0*SUP_HEXA/(DELTAZ/D[iz+1,nreg,grupo] + DELTAZ/D[iz,nreg,grupo]);
        diag[iz,nreg,grupo] := diag[iz,nreg,grupo] + EX;
        coeff[grupo,iz,nreg,8] := EX;
        end

      else if not ConMatrizDeRespuesta then begin
        EX := (SUP_HEXA*ALFA/(1.0 + DELTAZ*ALFA*0.5/D[iz,nreg,grupo]));
        coeff[grupo,iz,nreg,8] := EX;
        diag[iz,nreg,grupo] := diag[iz,nreg,grupo] + EX;
        end;

      end;

    diag0 := diag;

    for iz:=1 to NZ do
    for nreg:=1 to NTotalRegiones do
    for grupo:=1 to ng do
    if D[iz,nreg,grupo] > 0.0 then begin
      diag[iz,nreg,grupo] := diag0[iz,nreg,grupo] + SIGTOT[iz,nreg,grupo];
      end;

    for iz:=1 to NZ do
    for nreg:=1 to NTotalRegiones do
    if D[iz,nreg,1] > 0.0 then
    if TipoCalculo = directo then begin

      for grupo:=1 to ng do begin
        NUSIGF0[iz,nreg,grupo] :=NUSIGF0[iz,nreg,grupo]/KE0;
        NUSIGF[iz,nreg,grupo] :=NUSIGF[iz,nreg,grupo]/KE0;
        end;

      for grupo:=1 to ng do begin
        LS[grupo] := LS0[grupo]*FI[iz,nreg,grupo];

        for prec:=1 to NPREC do
          LS[grupo] := LS[grupo]
              + chis[prec,grupo]*Lambda[prec]/(1.0 + lambda[prec]*DT)*ConcPrecursor[prec,iz,nreg];

        LeftSide[iz,nreg,grupo] := (LS[grupo] + fuente)*VOL;

//         NUSIGF[iz,i,j,grupo] :=
//             NUSIGF[iz,i,j,grupo]*((1.0-BetaEfTotal) + BetaLambdaDT[grupo]);
         end


       end

     else if TipoCalculo = adiabatico then
       for grupo:=1 to ng do begin
         NUSIGF0[iz,nreg,grupo] :=NUSIGF0[iz,nreg,grupo]/KE0;
         NUSIGF[iz,nreg,grupo] :=NUSIGF[iz,nreg,grupo]/KE0;
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
 const
    XS: array[1..40] of real =
 (2.61698E+00, 1.36062E+00, 8.22178E-01, 4.40172E-01, 2.52215E-01,
  2.65897E-03, 6.43050E-02, 5.30847E-04, 1.24768E-11, 2.05254E-13,
  0.00000E+00, 1.43602E-03, 7.07928E-02, 4.15708E-06, 4.09713E-07,
  0.00000E+00, 0.00000E+00, 1.78107E-02, 5.22929E-02, 4.68823E-03,
  0.00000E+00, 0.00000E+00, 2.01360E-03, 5.36672E-02, 2.50017E-01,
  0.00000E+00, 0.00000E+00, 4.83469E-06, 3.51147E-01, 1.01826E-01,
  6.07205E-03, 5.86657E-04, 8.35188E-03, 6.70160E-02, 1.45868E-01,
  2.40071E-03, 2.52309E-04, 3.61949E-03, 2.90310E-02, 6.31894E-02);

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
  procedure Imprimir_XS_Red (NC1,NC2,NZ1,NZ2: integer);
  var
    nreg,kz,grupo,grupo1: integer;
    ss: real;

  begin
    for kz:=NZ1 to NZ2 do begin
      writeln(F);
      writeln(F,'XS DE LA RED NIVEL ',kz);
      for nreg:=NC1 to NC2 do begin
        writeln(F);
        writeln(F,'region ',nreg);

        for grupo:=1 to ng do begin
           if D[kz,nreg,grupo] = 0.0 then
             ss := 0.0
           else
             ss := 1.0/(3.0*D[kz,nreg,grupo]);

           write(F,FloatToStrF(ss,ffexponent,5,2),' ');
           end;

        writeln(F);
        for grupo:=1 to ng do begin
          for grupo1:=1 to ng do
            if grupo1 = grupo then
              write(F,FloatToStrF(SIGTOT[kz,nreg,grupo],ffexponent,5,2),' ')
            else
              write(F,FloatToStrF(SigSCATT[kz,nreg,grupo,grupo1],ffexponent,5,2),' ');
          writeln(F);
          end;

        for grupo:=1 to ng do
          write(F,FloatToStrF(NUSIGF[kz,nreg,grupo],ffexponent,5,2),' ');

        writeln(F);
        for grupo:=1 to ng do
          write(F,FloatToStrF(SIGF[kz,nreg,grupo],ffexponent,5,2),' ');

        end;
      end;

   Close(F);
   ShowMessage ('Listo');
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

      if num_dedos = 0 then
        if tipo = '1.8%' then NT := 1
        else if tipo = '3.1%' then
          if VQ = ' 0  ' then NT := 2
          else if VQ = ' 6VQ' then NT := 3
          else if VQ = '12VQ' then NT := 4;

      if num_dedos = 12 then
        if tipo = '1.8%' then NT := 5
        else if tipo = '3.1%' then
          if VQ = ' 0  ' then NT := 6
          else if VQ = ' 6VQ' then NT := 8 ;

      if num_dedos = 18 then
        if tipo = '3.1%' then
          if VQ = ' 0  ' then NT := 7
          else if VQ = ' 6VQ' then NT := 9 ;

      if (trozo = 6) or (trozo = 11) then NT := NT + 9 ;

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

  materiales ( 4,'3.1%','12VQ', 0, no_barra) ;
  materiales ( 5,'3.1%',' 6VQ',12, barra_E2) ;
  materiales ( 6,'3.1%','12VQ', 0, no_barra) ;

  materiales ( 7,'3.1%',' 0  ', 0, no_barra) ;
  materiales ( 8,'3.1%',' 0  ',18, barra_D3) ;
  materiales ( 9,'3.1%',' 6VQ',18, barra_F2) ;
  materiales (10,'3.1%',' 0  ', 0, no_barra) ;

  materiales (11,'3.1%',' 6VQ', 0, no_barra) ;
  materiales (12,'3.1%',' 0  ',18, barra_C4) ;
  materiales (13,'3.1%','12VQ', 0, no_barra) ;
  materiales (14,'3.1%',' 0  ',18, barra_G2) ;
  materiales (15,'3.1%',' 6VQ', 0, no_barra) ;

  materiales (16,'3.1%',' 0  ',18, barra_B5) ;
  materiales (17,'3.1%',' 6VQ',18, barra_D4) ;
  materiales (18,'3.1%',' 6VQ',0 , no_barra) ;
  materiales (19,'3.1%',' 0  ',18, barra_H2) ;

  materiales (20,'3.1%',' 0  ', 0, no_barra) ;
  materiales (21,'3.1%','12VQ', 0, no_barra) ;
  materiales (22,'3.1%',' 6VQ',12, barra_E4) ;
  materiales (23,'3.1%','12VQ', 0, no_barra) ;
  materiales (24,'3.1%',' 0  ', 0, no_barra) ;

  materiales (25,'3.1%',' 0  ',18, barra_B6) ;
  materiales (26,'3.1%',' 6VQ',18, barra_D5) ;
  materiales (27,'3.1%',' 6VQ',18, barra_F4) ;
  materiales (28,'3.1%',' 0  ',18, barra_H3) ;

  materiales (29,'3.1%','12VQ', 0, no_barra) ;
  materiales (30,'3.1%',' 6VQ', 0, no_barra) ;
  materiales (31,'1.8%',' 0  ',12, barra_E5) ;
  materiales (32,'3.1%',' 6VQ',18, barra_G4) ;
  materiales (33,'3.1%','12VQ', 0, no_barra) ;

  materiales (34,'3.1%',' 6VQ',18, barra_B7) ;
  materiales (35,'3.1%',' 6VQ',12, barra_D6) ;
  materiales (36,'3.1%',' 6VQ',12, barra_F5) ;
  materiales (37,'3.1%',' 0  ',18, barra_H4) ;

  materiales (38,'3.1%',' 6VQ', 0, no_barra) ;
  materiales (39,'3.1%','12VQ', 0, no_barra) ;
  materiales (40,'3.1%',' 6VQ',18, barra_E6) ;
  materiales (41,'3.1%','12VQ', 0, no_barra) ;
  materiales (42,'3.1%',' 6VQ', 0, no_barra) ;

  materiales (43,'3.1%',' 6VQ',12, barra_B8) ;
  materiales (44,'3.1%',' 6VQ',18, barra_D7) ;
  materiales (45,'3.1%',' 6VQ', 0, no_barra) ;
  materiales (46,'3.1%',' 6VQ',12, barra_H5) ;

  materiales (47,'3.1%',' 6VQ', 0, no_barra) ;
  materiales (48,'3.1%',' 0  ',18, barra_C8) ;
  materiales (49,'3.1%','12VQ', 0, no_barra) ;
  materiales (50,'3.1%',' 6VQ',18, barra_G6) ;
  materiales (51,'3.1%',' 6VQ', 0, no_barra) ;

  materiales (52,'3.1%',' 6VQ', 0, no_barra) ;
  materiales (53,'3.1%',' 0  ',18, barra_D8) ;
  materiales (54,'3.1%',' 0  ',18, barra_F7) ;
  materiales (55,'3.1%',' 6VQ', 0, no_barra) ;

  materiales (56,'3.1%','12VQ', 0, no_barra) ;
  materiales (57,'3.1%',' 0  ',18, barra_E8) ;
  materiales (58,'3.1%','12VQ', 0, no_barra) ;

  materiales (59,'3.1%',' 0  ', 0, no_barra) ;
  materiales (60,'3.1%',' 0  ', 0, no_barra) ;

  materiales (61,'3.1%',' 6VQ', 0, no_barra) ;

  end;


(*******************************************************************************)
  procedure InterpolateXS ;
  const
    NUMERO: array[1..18] of string[2] =
     ('01','02','03','04','05','06','07','08','09','10',
      '11','12','13','14','15','16','17','18');


   DBORO: array[1..40] of real =
  (-1.956000E-07, 3.594000E-07, 8.848000E-07, 1.182900E-06, 1.420900E-05,
    1.000000E-08,-2.600000E-07,-2.400000E-09, 4.000000E-17, 8.000000E-19,
    0           , 1.500000E-08, 7.400000E-07, 4.400000E-11, 5.100000E-12,
    0           , 0           , 8.200000E-07, 2.700000E-07, 2.500000E-08,
    0           , 0           , 2.130000E-07, 7.670000E-06,-6.100000E-06,
    0           , 0           , 1.000000E-12, 4.000000E-07, 1.741000E-05,
    2.800000E-08, 1.100000E-09, 3.100000E-08,-7.400000E-07,-3.000000E-07,
    9.000000E-09, 5.000000E-10, 1.300000E-08,-3.200000E-07,-1.100000E-07);



   DXENON: array[1..40] of real =
   (-2.547448E-20, 4.379789E-20, 2.495397E-20, 1.292825E-19, 3.780344E-18,
     9.408420E-22,-3.364018E-20,-3.045891E-22, 2.933532E-28, 4.888317E-30,
     0.000000E+00, 7.039394E-22, 9.036145E-20, 5.441993E-24, 6.328686E-25,
     0.000000E+00, 0.000000E+00, 8.054691E-21, 9.076756E-20, 7.668878E-21,
     0.000000E+00, 0.000000E+00, 2.542304E-20, 8.210370E-19,-6.843103E-19,
     0.000000E+00, 0.000000E+00, 5.658589E-24, 3.932584E-19, 2.137810E-18,
     3.519697E-21, 1.346961E-22, 6.085014E-21,-2.260051E-19,-6.152701E-19,
     1.150670E-21, 6.024096E-23, 2.633004E-21,-9.787464E-20,-2.662786E-19);

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
procedure agregar_xenon;
 var
   grupo, grupo1: integer;

 begin
   DELTA := ConcXe[canal,trozo] - ConcXe0;
   pos := 1;

   for grupo:=1 to ng do begin
     DIF[canal,trozo,grupo] := DIF[canal,trozo,grupo] + DELTA*DXENON[pos];
     inc(pos);
     end;

   for grupo:=1 to ng do for grupo1:=1 to ng do begin
     if grupo = grupo1 then
       SIGMABS[canal,trozo,grupo] :=
         SIGMABS[canal,trozo,grupo] + DELTA*DXENON[pos]
    else
      SCATTERING[canal,trozo,grupo,grupo1] :=
         SCATTERING[canal,trozo,grupo,grupo1] + DELTA*DXENON[pos];
    inc(pos);
    end;

  for grupo:=1 to ng do begin
    NUSIGMAFIS[canal,trozo,grupo] := NUSIGMAFIS[canal,trozo,grupo] +
         DELTA*DXENON[pos];
    inc(pos);
    end;

  for grupo:=1 to ng do begin
    SIGMAFIS[canal,trozo,grupo] := SIGMAFIS[canal,trozo,grupo] +
         DELTA*DXENON[pos];
    inc(pos);
    end;

  end; (* agregar_xenon *)

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

(* --------------------------------------------------------------------- *)
  begin
    ttt0 := time();
    asignar_materiales;

    if ArchivoIntMultiple <> '' then begin

      if not InterpolacionMultipleInicializada then
      for i:=1 to 18 do begin
        SetFileName(ArchivoIntMultiple + NUMERO[i] + '.TXT',i,3);
        InterpolacionMultipleInicializada := true;
        end;

      for canal:=1 to NCANALES do
      for trozo:=1 to NTROZOS do begin
        P1 := BurnUp[canal,trozo];
        P2 := TComb[canal,trozo];
        P3 := TempRefr[canal,trozo];
        P4 := DensRefr[canal,trozo];
        P5 := 0.0; //ConcXe[canal,trozo];
        P6 := 0.0; //PPMBoro[canal,trozo]*1.0E-4;
        P7 := 0.0;

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
          DELTA := ConcXe[canal,trozo] - ConcXe0;
          N00 := 241;
          agregar_Xenon (N00, DELTA);
          end;

        end;
   
     end

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

      if ConXenon then begin
        DELTA := ConcXe[canal,trozo] - ConcXe0;
        N00 := 241;
        agregar_Xenon (N00, DELTA);
        end;

      if ConSamario then begin
        DELTA := ConcSm[canal,trozo] - ConcSamario0;
        N00 := 281;
        agregar_Samario (N00, DELTA);
        end;

      if ConReacoplamientoTermohidraulico then begin
        N00 := 41;
        DELTA := TComb[canal,trozo] - TComb0;
        agregar (N00, DELTA);

        N00 := 81;
        DELTA := TempRefr[canal,trozo] - TempRefr0;
        agregar (N00, DELTA);


        if DensRefr[canal,trozo] <= 550.0 then begin
         Fact1 := (550.0 - DensRefr[canal,trozo])/150.0;
         Fact2 := (DensRefr[canal,trozo] - 400.0)/150.0;
         agregar (121,Fact1);
         agregar (161,Fact2);
         end

        else if DensRefr[canal,trozo] <= 675.0 then begin
          Fact1 := (675.0 - DensRefr[canal,trozo])/125.0;
          agregar (161,Fact1);
          end

        else begin
          Fact1 := (DensRefr[canal,trozo] - 675.0)/78.0;
          agregar (201,Fact1);
          end;

        end;

      if ConBoro then agregar_boro;
//        N00 := 161;
//        DELTA := PPMBoro[canal,trozo];
//        end;

      end;

    for canal:=1 to NCANALES do
    for trozo:=1 to NTROZOS do
      for grupo:=1 to ng do
        if DIF[canal,trozo, grupo] >= 0.0 then
          DIF[canal,trozo, grupo] := 1.0/(3.0*DIF[canal,trozo, grupo]);

  NUSIGMAFIS0 :=NUSIGMAFIS;

  if not ArchivoSalidaAbierto then begin
    assign(F,ArchivoDeSalida);
    rewrite(F);
    ArchivoSalidaAbierto := true;
    end;
{
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
 }
//   Close(F);
//   ShowMessage ('Listo');
//   Halt(0);

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

        if iz < NREFL1 then begin
          pos := 1;
          for grupo:=1 to ng do begin
            D[iz,canal,grupo] := XS_REFL_Inf[pos];
            inc(pos);
            end;

          for grupo:=1 to ng do for grupo1:=1 to ng do begin
            if grupo = grupo1 then
              SIGTOT[iz,canal,grupo] := XS_REFL_Inf[pos]
            else
              SIGSCATT[iz,canal,grupo,grupo1] :=  XS_REFL_Inf[pos];
            inc(pos);
            end;

          for grupo:=1 to ng do for grupo1:=1 to ng do
            SIGTOT[iz,canal,grupo] := SIGTOT[iz,canal,grupo]
                + SIGSCATT[iz,canal,grupo,grupo1];

          for grupo:=1 to ng do begin
            NUSIGF[iz,canal,grupo] :=  0.0;
            NUSIGF0[iz,canal,grupo] :=  0.0;
            SIGF[iz,canal,grupo] :=  0.0;
            end;

          end

        else if iz <= NREFL2 then begin
          trozo := (iz-NREFL1+1);
          D[iz,canal] := DIF[canal,trozo];
          for grupo:=1 to ng do begin
            SIGTOT[iz,canal,grupo] :=  SIGMABS[canal,trozo,grupo];
            for grupo1:=1 to ng do
              SIGTOT[iz,canal,grupo] := SIGTOT[iz,canal,grupo]
                + SCATTERING[canal,trozo,grupo,grupo1];

            end;

          SIGSCATT[iz,canal] :=  SCATTERING[canal,trozo];

          NUSIGF[iz,canal] :=  NUSIGMAFIS[canal,trozo];
          NUSIGF0[iz,canal] :=  NUSIGMAFIS0[canal,trozo];
          SIGF[iz,canal] :=  SIGMAFIS[canal,trozo];
          end

        else begin
          pos := 1;
          for grupo:=1 to ng do begin
            D[iz,canal,grupo] := XS_REFL_Sup[pos];
            inc(pos);
            end;
          for grupo:=1 to ng do for grupo1:=1 to ng do begin
            if grupo = grupo1 then
              SIGTOT[iz,canal,grupo] := XS_REFL_Sup[pos]
            else
              SIGSCATT[iz,canal,grupo,grupo1] :=  XS_REFL_Sup[pos];
            inc(pos);
            end;

          for grupo:=1 to ng do for grupo1:=1 to ng do
            SIGTOT[iz,canal,grupo] := SIGTOT[iz,canal,grupo]
                + SIGSCATT[iz,canal,grupo,grupo1];

          for grupo:=1 to ng do begin
            NUSIGF[iz,canal,grupo] :=  0.0;
            NUSIGF0[iz,canal,grupo] :=  0.0;
            SIGF[iz,canal,grupo] :=  0.0;
            end;

          end

        end;

      end;

    for iz:=1 to NZ do begin

      for canal:=62 to NTotalRegiones do begin

        if iz < NREFL1 then begin
          pos := 1;
          for grupo:=1 to ng do begin
            D[iz,canal,grupo] := XS_REFL_Rad[pos];
            inc(pos);
            end;
          for grupo:=1 to ng do for grupo1:=1 to ng do begin
            if grupo = grupo1 then
              SIGTOT[iz,canal,grupo] := XS_REFL_Rad[pos]
            else
              SIGSCATT[iz,canal,grupo,grupo1] :=  XS_REFL_Rad[pos];
            inc(pos);
            end;

          for grupo:=1 to ng do for grupo1:=1 to ng do
            SIGTOT[iz,canal,grupo] := SIGTOT[iz,canal,grupo]
                + SIGSCATT[iz,canal,grupo,grupo1];

          end

        else if iz > NREFL2 then begin
          pos := 1;
          for grupo:=1 to ng do begin
            D[iz,canal,grupo] := XS_REFL_Rad[pos];
            inc(pos);
            end;

          for grupo:=1 to ng do for grupo1:=1 to ng do begin
            if grupo = grupo1 then
              SIGTOT[iz,canal,grupo] := XS_REFL_Rad[pos]
            else
              SIGSCATT[iz,canal,grupo,grupo1] :=  XS_REFL_Rad[pos];
            inc(pos);
            end;

          for grupo:=1 to ng do for grupo1:=1 to ng do
            SIGTOT[iz,canal,grupo] := SIGTOT[iz,canal,grupo]
                  + SIGSCATT[iz,canal,grupo,grupo1];

          end

        else begin
          pos := 1;
          for grupo:=1 to ng do begin
            D[iz,canal,grupo] := XS_REFL_Rad[pos];
            inc(pos);
            end;

          for grupo:=1 to ng do for grupo1:=1 to ng do begin
            if grupo = grupo1 then
              SIGTOT[iz,canal,grupo] := XS_REFL_Rad[pos]
            else
              SIGSCATT[iz,canal,grupo,grupo1] :=  XS_REFL_Rad[pos];
            inc(pos);
            end;

          for grupo:=1 to ng do for grupo1:=1 to ng do
            SIGTOT[iz,canal,grupo] := SIGTOT[iz,canal,grupo]
                + SIGSCATT[iz,canal,grupo,grupo1];
          end;

        for grupo:=1 to ng do begin
          NUSIGF[iz,canal,grupo] :=  0.0;
          NUSIGF0[iz,canal,grupo] :=  0.0;
          SIGF[iz,canal,grupo] :=  0.0;
          end;

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
   Intercalado: boolean;

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

   while (not eof(Entrada)) do begin
     readln(Entrada,S);

     uppercase(S);

     if (S = '') or (S[1] = '!') then

     else if S[1] = '*' then
       if ArchivoDeSalida <> '' then
         writeln(F,S)
       else

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
        readln(Entrada,ArchivoIntMultiple)

     else if E('POTENCIA') and E('FISION') then
       readln (Entrada,PPorFision)

     else if E('TABLA') and E('VALORES') and E('CENTRALES') then
       readln (Entrada,ConcXE0, TComb0, DensRefr0, TempRefr0)

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
       for j:=1 to 10 do read (Entrada, XS_REFL_Inf[j]) ;
       for j:=1 to 10 do read (Entrada, XS_REFL_Sup[j]) ;
       for j:=1 to 10 do read (Entrada, XS_REFL_Rad[j]) ;
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
           read (Entrada,TCOMB[canal,trozo]);
         readln(Entrada);
         end;

       end

     else if E('TEMPERATURA') and E('REFRIGERANTE') then begin
       for canal:=1 to NCANALES do begin
         if intercalado then
           readln(Entrada);

         for trozo:=1 to NTROZOS do
           read (Entrada,TempRefr[canal,trozo]);
         readln(Entrada);
         end;
       end

     else if E('DENSIDAD') and E('REFRIGERANTE') then begin
       for canal:=1 to NCANALES do begin
         if intercalado then
           readln(Entrada);

         for trozo:=1 to NTROZOS do
           read (Entrada,DensRefr[canal,trozo]);
         readln(Entrada);
         end;
       end

     else if E('QUEMADO') then
       for canal:=1 to NCANALES do begin
         if intercalado then
           readln(Entrada);

         for trozo:=1 to NTROZOS do
           read (Entrada,BurnUp[canal,trozo]);
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
   canal,trozo, prec,grupo, nreg, iz, i, i1, i2, j1, j2: integer;

 begin
//  ArchivoSalidaAbierto := false;
  ConSamario := false;
  ConXenon := false;                   ConReacoplamientoTermohidraulico := false;
  ImprimirIteraciones := false;        ImprimirNumIteraciones := false;
  ConBarras := false;                  ConQuemado := false;
  ParaCompararConPUMA := false;        ImprimirIteraciones := false;
  ImprimirNumIteraciones := false;     ImprimirDistribucion := false;
  ImprimirPorCanal := false;
  UNIDAD := 1.0;                       FactPuma := 1.0;
  PrecisionEstacionario := 1.0E-6;     PrecisionCinetico := 1.0E-6;
  ConcXE0 := 1.3056E+15;               ConcSamario0 := 4.4004E+16; // 1.0856E16;
  TComb0 := 450.0;
  DensRefr0 := 675.0;                  TempRefr0 := 305.0;
  PuedeImprimir := true;               PPorFision := 1.0;
  MaxIteraciones := 600;               variacion_lineal := false;

  ArchivoIntMultiple := '';            InterpolacionMultipleInicializada := false;
  ConMatrizDeRespuesta := false;
  BETA0 := 1.6;                        interpolacion_lineal := false;

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
  for prec:=1 to nprec do
    BetaNuclearTotal := BetaNuclearTotal + BetaNuclear[prec];

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
    end;

  for canal:=1 to NCANALES do
    for trozo:=1 to NTROZOS do
      NumMaterial[canal,trozo] := 2;

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
  for nreg:=1 to NTotalRegiones do
  for grupo:=1 to ng do begin
    FI[iz,nreg,grupo] := 1.0;
    FIA[iz,nreg,grupo] := 1.0;
    end;

  VOLCeldilla := 6.0*DELTAZ*ALTURA*ALTURA/sqrt(3);
  VolTrozo := VOLCeldilla;

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
   iz,grupo: integer;
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
       for grupo:=1 to ng do
         flujo[canal,trozo,grupo] := flujo[canal,trozo,grupo]
              + FI[iz,canal,grupo]*VolCeldilla;

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
   canal,trozo,grupo,sxt,iz,nreg: integer;
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
  for nreg:=1 to NTotalRegiones do
  for grupo:=1 to ng do
    FI[iz,nreg,grupo] := FI[iz,nreg,grupo]*FNORM;

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
   FACT: real;
//   vv: array[1..91] of integer;
//   i: integer;

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
//    writeln (F,'    Temp Comb Media = ',TCombMedia:1:1);

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
      PotTotal/(NCANALES*NTROZOS)*FForma:1:3,
          '   en el canal: ',CanalPMax:3,'  trozo: ',TrozoPMax:3);

  if DT = 0.0 then begin
    writeln (F,'     KE (estático) = ',KE0:1:5);
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

    writeln(F,'Potencias por canal');

    for canal:=1 to NCANALES do begin
      write (F,'CAN',canal:3,' ',potcan[canal]:1:4,'  ');
      if canal mod 6 = 0 then writeln(F);
      end;

   end;

 if ParaCompararConPUMA then begin
   writeln(F);
   writeln(F);

   FACT := (PotTotal-PotResidual)/PotTotal;

   for canal:=1 to NCANALES do begin
     write (F,' ',potcan[canal]*FACT:8:4,' ');
     if canal mod 8 = 0 then writeln(F);
     end;

   end;

 if ImprimirDistribucion then begin

   writeln(F);
   writeln(F);

   for canal:=1 to NCANALES do begin
     writeln(F);
     writeln (F,'canal ',canal,'   ',potcan[canal]:1:4);
     for trozo:=1 to NTROZOS do begin
        write(F,POT[canal,trozo]/VOLTROZO*1.0E6:8:2,' ');
        if trozo mod 10 = 0 then writeln(F);
        end;
     end;
 (*
   for canal:=1 to NCANALES do begin
     writeln(F);
     writeln(F);
     writeln (F,'canal ',canal);
     for trozo:=1 to NTROZOS do begin
        write(F,BurnUp[canal,trozo]:8:2,' ');
        if trozo mod 10 = 0 then writeln(F);
        end;
     end;


   for canal:=1 to NCANALES do begin
     writeln(F);
     writeln(F);
     writeln (F,'canal ',canal);
     for trozo:=1 to NTROZOS do begin
        write(F,NumMaterial[canal,trozo]:8,' ');
        if trozo mod 10 = 0 then writeln(F);
        end;
     end;
   *)
    end;

  writeln(F);
(*
  for canal:=1 to 91 do
      vv[canal] := 0;

  for canal:=1 to 91 do
    for i:=1 to 6 do
    if vecino[canal,i] <> 0 then
      vv[vecino[canal,i]] := vv[vecino[canal,i]] + 1;

  writeln(F);
  for canal:=1 to 91 do
    writeln (F,'canal ',canal,'  ',vv[canal]);

    writeln(F);

  for canal:=1 to 91 do begin
    writeln (F,'canal ',canal);
    for i:=1 to 6 do
    if vecino[canal,i] <> 0 then
      writeln (F,vecino[vecino[canal,i],1+(i+2) mod 6],'  ',canal);
    end;

  *)
  end;

//**************************************************************************/
//*  Resuelve el problema de autovalores de flujo y K efectivo para la     */
//*  matriz hallada en estado estacionario                                 */
//**************************************************************************/

 procedure Iteraciones (TipoCalculo: TipoCalculoType; precision: real);
 var
   nreg,iz, iteracion, NNN, n, NP, MaxIter, grupo, grupo1, nborde, nvecino: integer;
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
  for nreg:=1 to NTotalRegiones do
  for iz:=1 to NZ do
  if D[iz,nreg,1] > 0.0 then
    for grupo:=1 to ng do
      ACC0 := ACC0 + FI[iz,nreg,grupo]

  else
    for grupo:=1 to ng do begin
      FI[iz,nreg,grupo] := 0.0;
      FIA[iz,nreg,grupo] := 0.0;
      end;

  OMEGA0 := 1.0;
  OMEGA := 1.0;
  NNN := 0;

  while (iteracion <= MaxIter) and (EPS >= Precision) do begin
    ACC := 0.0;

    for iz:=1 to NZ do
    for nreg:=1 to NTotalRegiones do begin

      for grupo:=1 to ng do
        FLU[grupo] := FI[iz,nreg,grupo];

      for grupo:=1 to ng do begin
        TermIzquierdo := 0.0;

        if D[iz,nreg,1] > 0.0 then
          Diagonal := diag[iz,nreg,grupo]

        else
          Diagonal := 0.0;

        if D[iz,nreg,1] > 0.0 then begin
          for nborde := 1 to 6 do begin
            nvecino := vecino[nreg,nborde];

            if (nvecino > 0) and (D[iz,nvecino,1] > 0.0) then
              TermIzquierdo := TermIzquierdo + coeff[grupo,iz,nreg,nborde]*FI[iz,nvecino,grupo];

            end;

          if (iz > 1) and (D[iz-1,nreg,1] > 0.0) then
            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,nreg,7]*FI[iz-1,nreg,grupo];

          if (iz < NZ) and (D[iz+1,nreg,1] > 0.0) then
            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,nreg,8]*FI[iz+1,nreg,grupo];

          for grupo1:=1 to ng do
            TermIzquierdo := TermIzquierdo + espectro_estatico[grupo]*NUSIGF[iz,nreg,grupo1]*FLU[grupo1]/KE;

          for grupo1:=1 to ng do
            if grupo1 <> grupo then
            TermIzquierdo := TermIzquierdo + SIGSCATT[iz,nreg,grupo1,grupo]*FLU[grupo1];

          TermIzquierdo := TermIzquierdo/Diagonal;

          FLU[grupo] := FI[iz,nreg,grupo] + BETA0*(TermIzquierdo - FI[iz,nreg,grupo]);
          end;

        end;

      for grupo:=1 to ng do
        ACC := ACC + FLU[grupo];

      for grupo:=1 to ng do
        FI[iz,nreg,grupo] := FLU[grupo];

      end;

    KE := KE*ACC/ACC0;
    EPSKE := abs(KE-KE_ANT)/KE;
    KE_ANT := KE;
    ACC0 := ACC;
    EPS := 0.0;

    NP := 0;

    for nreg:=1 to NTotalRegiones do
    for iz:=NZMIN to NZMAX do
    if D[iz,nreg,1] > 0.0 then
      for grupo:=1 to ng do begin
      ex := abs((FI[iz,nreg,grupo] - FIA[iz,nreg,grupo])/FI[iz,nreg,grupo]);
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
      for nreg:=1 to NTotalRegiones do
      for iz:=NZMIN to NZMAX do
      for grupo:=1 to ng do begin
        FI[iz,nreg,grupo] := FIA[iz,nreg,grupo] + OMEGA*(FI[iz,nreg,grupo] - FIA[iz,nreg,grupo]);
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

//**************************************************************************/
//*  Resuelve el problema de cálculo del flujo para la matriz hallada      */
//*  en estado no estacionario para problemas cinético espaciales          */
//**************************************************************************/

 procedure CalculoCineticoEspacial(precision: real);
 var
   iteracion, nreg, NNN, n, i,j,iz, NP, grupo, grupo1, nborde, nvecino: integer;
   MaxIter: integer;
   EPS0: real;
   Diagonal, TermIzquierdo: real;
   FLU: array[1..ng] of real;
   OMEGA, EPSOMEGA, OMEGA0, OMEGAX: real;
   KE, EPS, ALFA, ACCF, ACCA, VOL, ex: real;
   TTT00: TDateTime;


 begin
  TTT00 := Time;
  VOL := VolCeldilla;

  if (ImprimirIteraciones or ImprimirNumIteraciones)
                 and not ArchivoSalidaAbierto then begin
    assign(F,ArchivoDeSalida);
    rewrite(F);
    ArchivoSalidaAbierto := true;
    end;

  if ImprimirIteraciones and PuedeImprimir then begin
    Writeln(F);
    Writeln(F,'Cálculo cinético espacial de flujo y K efectivo');
    Writeln(F,'BETA = ',BETA0:5:2);
    Writeln(F);
    Writeln(F,' IT  Conv Flujo  Kefectivo  Sigma    Omega');
    end;

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

    for nreg:=1 to NTotalRegiones do
    for iz:=1 to NZ do begin

      for grupo:=1 to ng do
        FLU[grupo] := FI[iz,nreg,grupo];

      for grupo:=1 to ng do begin
        TermIzquierdo := 0.0;

        if D[iz,nreg,1] > 0.0 then
          Diagonal := diag[iz,nreg,grupo]

        else begin
          Diagonal := 0.0;
          FI[iz,nreg,grupo] := 0.0;
          FIA[iz,nreg,grupo] := 0.0;
          end;

        if D[iz,nreg,1] > 0.0 then begin

          for nborde := 1 to 6 do begin
            nvecino := vecino[nreg,nborde];

            if (nvecino > 0) and (D[iz,nvecino,1] > 0.0) then
              TermIzquierdo := TermIzquierdo + coeff[grupo,iz,nreg,nborde]*FI[iz,nvecino,grupo];

            end;

          if (iz > 1) and (D[iz-1,nreg,1] > 0.0) then
            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,nreg,7]*FI[iz-1,nreg,grupo];

          if (iz < NZ) and (D[iz+1,nreg,1] > 0.0) then
            TermIzquierdo := TermIzquierdo + coeff[grupo,iz,nreg,8]*FI[iz+1,nreg,grupo];

          ex := TermIzquierdo;

          for grupo1:=1 to ng do
            TermIzquierdo := TermIzquierdo + espectro_dinamico[grupo]*NUSIGF[iz,nreg,grupo1]*FLU[grupo1];

          for grupo1:=1 to ng do
            if grupo1 <> grupo then
              TermIzquierdo := TermIzquierdo + SIGSCATT[iz,nreg,grupo1,grupo]*FLU[grupo1];

          FLU[grupo] := (TermIzquierdo + LeftSide[iz,nreg,grupo])/(Diagonal + LS0[grupo]*VOL);

          FLU[grupo] := FIA[iz,nreg,grupo] + BETA0*(FLU[grupo] - FIA[iz,nreg,grupo]);

          for grupo1:=1 to ng do
            ACCF := espectro_estatico[grupo]*NUSIGF0[iz,nreg,grupo1]*FLU[grupo1]*VOL + ACCF;

          ACCA := - ex + diagonal*FLU[grupo] + ACCA;

          for grupo1:=1 to ng do
            ACCA := - SIGSCATT[iz,nreg,grupo1,grupo]*FLU[grupo1]*ord(grupo1 <> grupo)  + ACCA;
          end;  //  D

        end;  // grupo

      for grupo:=1 to ng do
        FI[iz,nreg,grupo] := FLU[grupo];

      end; //iz

    KE := ACCF/ACCA;
    EPS := 0.0;
    NP := 0;

    for iz:=1 to NZ do
    for nreg:=1 to NTotalRegiones do
    if D[iz,nreg,1] > 0.0 then begin
      inc(NP);
      for grupo:=1 to ng do
        EPS := EPS + abs((FI[iz,nreg,grupo] - FIA[iz,nreg,grupo])/FI[iz,nreg,grupo])
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
      for nreg:=1 to NTotalRegiones do
      for iz:=1 to NZ do
      for grupo:=1 to ng do begin
        FI[iz,nreg,grupo] := FIA[iz,nreg,grupo] + OMEGA*(FI[iz,nreg,grupo] - FIA[iz,nreg,grupo]);
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
    writeln (F);
    writeln (F,'Método Directo:');
    writeln  (F,'Tiempo para iteraciones: ',TTT00:1:3);
    writeln (F,'Num Iteraciones: ', NumIteraciones);
    writeln (F);
    end;

  RO := (KE-1)/KE;

  end; (*  CalculoCineticoEspacial  *)

(*******************************************************************************)
 procedure precursores(DT: real);
 var
   iz,prec,grupo,nreg: integer;
   FF: real;

 begin
   for iz:=1 to NZ do
   for nreg:=1 to NTotalRegiones do
   if D[iz,nreg,1] > 0.0 then begin
     FF := 0.0;
     for grupo:=1 to ng do
       FF := FF + NUSIGF0[iz,nreg,grupo]*FIANT[iz,nreg,grupo];

     for prec:=1 to NPREC do
     if DT = 0.0 then
       ConcPrecursor[prec,iz,nreg] := FF*BetaLambda[prec]/KE0
     else
       ConcPrecursor[prec,iz,nreg] :=
         (ConcPrecursor[prec,iz,nreg] + FF*BetaNuclear[prec]*DT)/(1 + Lambda[prec]*DT);
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
   prec,iz,nreg, grupo: integer;
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
   for nreg:=1 to NTotalRegiones do
   if D[iz,nreg,1] > 0.0 then begin
     for grupo:=1 to ng do begin
       ACC := ACC + FI[iz,nreg,grupo]/VEL[grupo];
       ACCF := ACCF + FI[iz,nreg,grupo]*NUSIGF0[iz,nreg,grupo]/KE0;
       end;

     for prec:=1 to NPREC do
       conc[prec] := conc[prec] + ConcPrecursor[prec,iz,nreg];

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

   if variacion_lineal then
     linear_variation(cn, conc, der_log, fuen, r0,r1,DT)
   else
     Backwards_differences (cn, conc, der_log, RO, 1.0E6, DT);

   DerivadaLogaritmica := der_log;

   ACC := 0.0;

   for iz:=1 to NZ do
   for nreg:=1 to NTotalRegiones do
   if D[iz,nreg,1] > 0.0 then
     for grupo:=1 to ng do
       ACC := ACC + FI[iz,nreg,grupo]/VEL[grupo];

   ACC := cn/ACC;

   for iz:=1 to NZ do
   for nreg:=1 to NTotalRegiones do
   if D[iz,nreg,1] > 0.0 then
     for grupo:=1 to ng do
       FI[iz,nreg,grupo] := FI[iz,nreg,grupo]*ACC;

   potencias(0.0,DT);
   precursores(DT);
   if ConXenon then XENON(DTXenon);
   if ConSamario then SAMARIO(DTXenon);
   end;

(*******************************************************************************)
 procedure DirectMethod(DT, DTXenon: real);
 var
   iz,nreg, grupo: integer;
   CONC1,CONC2: real;

 begin
   CONC1 := 0.0;

   for iz:=1 to NZ do
   for nreg:=1 to NTotalRegiones do
   if D[iz,nreg,1] > 0.0 then
     for grupo:=1 to ng do
       CONC1 := CONC1 + FI[iz,nreg,grupo]/VEL[grupo];

   FIANT := FI;
   InterpolateXS;
   MeshContants (directo,DT);
   CalculoCineticoEspacial(PrecisionCinetico);
   potencias(0.0,DT);
   precursores(DT);

   CONC2 := 0.0;

   for iz:=1 to NZ do
   for nreg:=1 to NTotalRegiones do
   if D[iz,nreg,1] > 0.0 then
     for grupo:=1 to ng do
       CONC2 := CONC2 + FI[iz,nreg,grupo]/VEL[grupo];

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
   MeshContants (estatico,0.0);
   iteraciones(estatico, PrecisionEstacionario);
   potencias(POWER,0.0);
   FIANT := FI;
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
