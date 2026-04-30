import sys
import json
import time
from datetime import datetime
from itertools import combinations
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import substitution_matrices
from Bio.Align import PairwiseAligner


class Pipeline:
    def __init__ (self):
        """
        Inicialización de distintos atributos para almacenar datos y resultados intermedios.
        """
        self.sequences = []
        self.metadata = {}
        self.config = {}
        self.results = {}
        self.performance_results = {}

    def obtener_hora(self):
        """
        Nos devuelve la hora actual de ejecución.
        """

        return datetime.now().strftime("%H:%M:%S")
    
    # Primer paso.
    def read_config(self, ruta="config.json"): # Si no recibo una ruta, se fijará en que sea config.json.
        """
        Recibe la ruta a un fichero JSON y almacena su contenido en el atributo config.

        ------------
        Params:
            ruta: Ruta del fichero json.

        Return:
            config: Diccionario que contiene los parámetros definidos en el fichero.
        ------------
        """

        print(f"[{self.obtener_hora()}] Paso 1: Recibe la ruta a un fichero JSON y almacena su contenido en el atributo config.\n")

        try: # Añado una excepción por si el fichero no existe.
            with open(ruta) as f:
                self.config = json.load(f)

                # Comprobaciones de que puedo acceder perfectamente a los parámetros del diccionario config. 
                print(f'Fichero de entrada, input_file: {self.config["input_file"]}') # Fichero de entrada: seqs_file.fasta
                print(f'Formato de entrada, input_format: {self.config["input_format"]}') # Formato de entrada: fasta
                print(f'Fichero de salida, output_file: {self.config["output_file"]}') # Fichero de salida: salida.fasta
                print(f'Formato de salida, output_format: {self.config["output_format"]}') # Formato de salida: fasta
                print(f'Clasificación de las secuencias: {self.config["classify_sequences"]}')
                print(f'Cálculo de estadísticas descriptivas: {self.config["compute_basic_stats"]}')
                print(f'Filtrado por tamaño: {self.config["filter_by_length"]}\n')
                
                print("Recepción y almacenamiento correcto.\n")

                return self.config
                
        except FileNotFoundError:
            print("El fichero de configuración no existe.\n")
            return 

        except json.JSONDecodeError:
            print("El archivo JSON está mal formado.\n")
            return 

    # Segundo paso.
    def load_sequences(self, fichero, input_format):
        """
        Carga secuencias desde un fichero utilizando SeqIO.parse y las almacena en el estado interno del pipeline (self.sequences) como objetos SeqRecord.

        Además, inicializa el diccionario self.metadata con información básica de las secuencias cargadas, incluyendo la longitud mínima, máxima y media, así como los 
        contadores de secuencias por tipo biológico (ADN, ARN y proteína), que se inicializan a cero para su posterior cálculo.

        ------------
        Params:
            fichero: Fichero de secuencias.
            input_format: Formato de entrada.
        ------------
        """

        print(f"[{self.obtener_hora()}] Paso 2: Carga y lectura de un fichero, junto con el almacenaje de las secuencias cargadas.\n")

        try: # Añado una excepción por si el archivo no existe.
            for seq_record in SeqIO.parse(fichero, input_format):
                self.sequences.append(seq_record)
        except FileNotFoundError:
            print("El fichero de secuencias no existe.\n")
            return

        # Incializar self.metadata
        # Inicializar longitud de secuencias, longitud mínima, máxima y media.
        longitudes = [len(n) for n in self.sequences]

        # Para que no se produzca un error sino hay longitudes.
        if longitudes == []:
            print("No se pueden calcular las estadísticas porque no hay secuencias válidas.\n")
            return
        
        self.metadata["num_seqs"] = len(self.sequences)
        self.metadata["min_length"] = min(longitudes)
        self.metadata["max_length"] = max(longitudes)
        self.metadata["mean_length"] = sum(longitudes)/len(longitudes)
    
        # Inicializar el número de secuencias de cada tipo biológico.
        self.metadata["n_dna"] = 0
        self.metadata["n_rna"] = 0
        self.metadata["n_prot"] = 0

        print("Carga, lectura y almacenamiento correctas.\n")

    # Tercer paso.
    def process(self):
        """
        Realiza un procesamiento básico sobre las secuencias almacenadas en el pipeline.

        Calcula la longitud de cada secuencia y la almacena en self.metadata bajo la clave "seq_length", permitiendo disponer de información estructurada para análisis 
        posteriores.
        """

        print(f"[{self.obtener_hora()}] Paso 3: Realizar operaciones de análisis simple sobre las secuencias.\n")

        longitudes = []
        for n in self.sequences:
            longitudes.append(len(n))

        self.metadata["seq_length"] = longitudes

        print("Realización de operaciones de análisis correctas.\n")

    # Cuarto paso.
    def filter_by_length(self):
        """
        Filtra las secuencias del pipeline eliminando aquellas cuya longitud es inferior a un umbral mínimo definido internamente.

        El método modifica directamente el estado interno (self.sequences), manteniendo únicamente las secuencias válidas. 

        Además, actualiza las estadísticas de longitud en self.metadata y reinicia los contadores de tipos biológicos para su posterior recalculado.
        """

        print(f"[{self.obtener_hora()}] Paso 4: Elimina las secuencias que no cumplen cierto valor de longitud.\n")
        
        min_length = 5
        nuevas_secuencias = []

        for n in self.sequences:
            if len(n) >= min_length:
                nuevas_secuencias.append(n)

            else:
                print(f"Se elimina la secuencia {n}.\n")

        self.sequences = nuevas_secuencias 

        # Actualización de los datos de self.metadata
        print("Paso 4.1: Actualización de los datos\n")
        longitudes = [len(n) for n in nuevas_secuencias]

        if longitudes == []:
            print("No se pueden calcular las estadísticas porque no hay secuencias válidas.\n")
            return
    
        self.metadata["min_length"] = min(longitudes)
        self.metadata["max_length"] = max(longitudes)
        self.metadata["mean_length"] = sum(longitudes)/len(longitudes)
    
        # Actaulizar el número de secuencias de cada tipo biológico, las inicializamos a 0 de nuevo.
        self.metadata["n_dna"] = 0
        self.metadata["n_rna"] = 0
        self.metadata["n_prot"] = 0

        print("Filtrado y actualización de los datos correcto.\n")

    # Quinto paso.
    def classify_and_normalize (self):
        """
        Clasifica cada secuencia según su tipo biológico (ADN, ARN o proteína) en función de los caracteres que contiene.

        Para ello, compara el conjunto de letras de cada secuencia con conjuntos predefinidos de caracteres válidos. 
        Las secuencias de ARN se convierten automáticamente a ADN mediante back_transcribe(), mientras que las de ADN y proteínas se mantienen sin modificar.

        Además, actualiza en self.metadata el número de secuencias de cada tipo biológico y descarta aquellas que contienen caracteres no válidos.
        """

        print(f"[{self.obtener_hora()}] Paso 5: Clasificación por el tipo biológico de cada secuencia.\n")

        n_ADN = {"A", "C", "G", "T"}
        n_ARN = {"A", "C", "G", "U"}
        n_proteina = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y" }

        nuevas_secuencias = []

        for n in self.sequences:
            letras = set(n.seq) # Obtengo todas las letras sin repetir.
            if letras <= n_ADN:
                self.metadata["n_dna"] +=1
                nuevas_secuencias.append(n)

            elif letras <= n_ARN:
                self.metadata["n_rna"] +=1
                adn = n.back_transcribe()
                nuevas_secuencias.append(adn)

            elif letras <= n_proteina:
                self.metadata["n_prot"] +=1
                nuevas_secuencias.append(n)

            else:
                print(f"La secuencia {n} es inválida.\n")

        self.sequences = nuevas_secuencias 

        print("Clasificación de las secuencias según su tipo biológico correcta.\n")

    # Sexto paso.
    def compute_basic_stats(self):
        """
        Calcula estadísticas descriptivas globales sobre las secuencias almacenadas en el pipeline.

        Incluye el número total de secuencias, así como la longitud mínima, máxima y media. Estos valores se muestran por pantalla y se almacenan en self.metadata.

        También muestra el número de secuencias de cada tipo biológico previamente calculado durante la fase de clasificación.
        """

        print(f"[{self.obtener_hora()}] Paso 6: Cálculo de las estadísticas descriptivas globales sobre las secuencias.\n")

        # Número total de secuencias.
        print(f"Número total de secuencias: {len(self.sequences)}.\n")

        # Longitud mínima, máxima y media.
        longitudes = [len(n) for n in self.sequences]

        # Para que no se produzca un error sino hay longitudes.
        if longitudes == []:
            print("No se pueden calcular las estadísticas porque no hay secuencias válidas.\n")
            return
        
        self.metadata["min_length"] = min(longitudes)
        self.metadata["max_length"] = max(longitudes)
        self.metadata["mean_length"] = sum(longitudes)/len(longitudes)

        print(f"Secuencia de menor tamaño: {min(longitudes)}.")
        print(f"Secuencia de mayor tamaño: {max(longitudes)}.")
        print(f"Tamaño medio de la longitud de las secuencias: {sum(longitudes)/len(longitudes)}.\n")
    
        # Número de secuencias de cada tipo biológico.
        print(f'Cantidad de secuencias de ADN: {self.metadata["n_dna"]}')
        print(f'Cantidad de secuencias de ARN: {self.metadata["n_rna"]}')
        print(f'Cantidad de secuencias de proteínas: {self.metadata["n_prot"]}.\n')

        print("Estadísticas calculadas.\n")

    # Séptimo paso.
    def smith_waterman(self, seq1, seq2, matrix_name, gap_penalty):
        """
        Implementa el algoritmo de alineamiento local Smith-Waterman para dos secuencias biológicas.

        El método construye una matriz de puntuación utilizando programación dinámica, aplicando una matriz de sustitución (PAM, BLOSUM, etc.) y una penalización por 
        gaps.

        Posteriormente, realiza un proceso de backtracking desde la celda de mayor puntuación para reconstruir el mejor alineamiento local.

        ------------
        Params:
            seq1: Secuencia 1.
            seq2: Secuencia 2.
            matrix_name: Nombre de la matriz de sustitución que se utiliza, en este caso PAM250.
            gap_penalty: Puntuación de penalización de gaps.

        Return:
            alineamiento_seq1: Secuencia 1 alineada.
            alineamiento_seq2: Secuencia 2 alineada.
            puntuacion: Puntuación máxima de alineamiento local.
            matchs: Coincidencia exacta en la alineación.
            mismatchs: No coincidencia en la alineación.
            gaps: Insercción  o deleción en la alineación
        ------------
        """

        # Cargo la matriz de sustitución PAM/BLOSUM.
        matriz_S = substitution_matrices.load(matrix_name) # Aquí se está cargando la matriz PAM250 definida en el archivo config.json.

        # Creo la matriz de puntuación.
        fila = len(seq1) + 1
        columna = len(seq2) + 1

        matriz_N = []
        for i in range(fila):
            nueva_fila = []
            for j in range(columna):
                nueva_fila.append(0)

            matriz_N.append(nueva_fila)

        # Creo la matriz del camino.
        fila = len(seq1) + 1
        columna = len(seq2) + 1

        matriz_T = []
        for i in range(fila):
            nueva_fila = []
            for j in range(columna):
                nueva_fila.append("0")

            matriz_T.append(nueva_fila) 
        
        # Relleno las matrices.
        for i in range(1, len(seq1) + 1):
            for j in range(1, len(seq2) + 1):
                diagonal = matriz_N[i-1][j-1] + matriz_S[seq1[i-1], seq2[j-1]]
                arriba = matriz_N[i-1][j] + gap_penalty
                izquierda = matriz_N[i][j-1] + gap_penalty

                matriz_N[i][j] = max(0, diagonal, arriba, izquierda)

                if matriz_N[i][j] == 0:
                    matriz_T[i][j] = "0"
                
                elif matriz_N[i][j] == diagonal:
                    matriz_T[i][j] = "DIAG"

                elif matriz_N[i][j] == arriba:
                    matriz_T[i][j] = "UP"    

                else:
                    matriz_T[i][j] = "LEFT"

        # Localizo el máximo score.
        puntuacion = 0
        max_i = 0
        max_j = 0

        for i in range(fila):
            for j in range(columna):
                if matriz_N[i][j] > puntuacion:
                    puntuacion = matriz_N[i][j]
                    max_i = i
                    max_j = j 

        # Realizo el backtracking desde el valor máximo hasta llegar al 0.
        alineamiento_seq1 = ""
        alineamiento_seq2 = ""

        # Inicializo a cero los contadores para la puntuación de "matchs", la puntuación de "mismatchs" y la puntuación de "gaps".
        matchs = 0
        mismatchs = 0
        gaps = 0

        i = max_i
        j = max_j


        while matriz_N[i][j] != 0:
            if matriz_T[i][j] == "DIAG":
                alineamiento_seq1 += seq1[i-1]
                alineamiento_seq2 += seq2[j-1]

                # Sumo en los contadores.
                if seq1[i-1] == seq2[j-1]:
                    matchs += 1
                else: 
                    mismatchs += 1

                i -= 1
                j -= 1
                
            elif matriz_T[i][j] == "UP":
                alineamiento_seq1 += seq1[i-1]
                alineamiento_seq2 += "-"
                gaps += 1
                i -= 1

            elif matriz_T[i][j] == "LEFT":
                alineamiento_seq1 += "-"
                alineamiento_seq2 += seq2[j-1]
                gaps += 1
                j -= 1
            
            else: 
                break

        # Invierto el alineamiento.
        alineamiento_seq1 = alineamiento_seq1[::-1]
        alineamiento_seq2 = alineamiento_seq2[::-1]

        return alineamiento_seq1, alineamiento_seq2, puntuacion, matchs, mismatchs, gaps

    def medir_tiempo_smith_waterman_propia(self):
        """
        Mide el tiempo de ejecución del algoritmo Smith-Waterman (nuestra implementación) sobre todas las combinaciones de pares de secuencias sin repetir.

        Genera todos los pares posibles de secuencias. Para cada par alinea y mide el tiempo que tarda e imprime el tiempo total y los tiempos de cada par.

        ------------
        Return:
            tiempos_pares: Diccionario que contiene los tiempos de cada par.
            tiempo_total: Tiempo total de ejecución.
        ------------
        """
 
        print("Medición de tiempo: Smith-Waterman propia.")

        # Diccionario donde guardamos el tiempo de cada par.
        tiempos_pares = {}
 
        #Verificamos que hay al menos 2 secuencias.
        if len(self.sequences) < 2:
            print("Error: Se necesitan al menos 2 secuencias para comparar.\n")
            return {}, 0
 
        # Obtenemos todas las combinaciones de pares de secuencias sin repetir.
        pares_secuencias = list(combinations(range(len(self.sequences)), 2))
        print(f"Número total de pares a alinear: {len(pares_secuencias)}.\n")
        
        # Medimos el tiempo de inicio.
        tiempo_inicio_total = time.time()
        
        #Ejecutamos el algoritmo de smith-waterman para cada par.
        for indice_seq1, indice_seq2 in pares_secuencias:
            seq1 = str(self.sequences[indice_seq1].seq)
            seq2 = str(self.sequences[indice_seq2].seq)

            # Tiempo inicio del par.
            inicio_par = time.time()

            self.smith_waterman(seq1, seq2, self.config["matrix_name"], self.config["gap_penalty"])

            # Tiempo fin del par.
            fin_par = time.time()

            # Tiempo total del par.
            tiempo_par_total = fin_par - inicio_par

            # Actulizamos el diccionario del tiempo de cada par.
            tiempos_pares[(indice_seq1, indice_seq2)] = tiempo_par_total

            print(f"Tiempo del par del algoritmo propio ({indice_seq1}, {indice_seq2}) es {tiempo_par_total:.4f} segundos.")

        # Registramos tiempo de fin y calculamos el total.
        tiempo_fin_total = time.time()
        tiempo_total = tiempo_fin_total - tiempo_inicio_total

        print(f"Tiempo total del algoritmo propio: {tiempo_total:.4f} segundos.\n")
 
        return tiempos_pares, tiempo_total
 
    def medir_tiempo_smith_waterman_biopython(self):
        """
        Mide el tiempo de ejecución del algoritmo Smith-Waterman (versión BioPython) sobre todas las combinaciones de pares de secuencias sin repetir.

        Genera todos los pares posibles de secuencias. Para cada par alinea con Biopython y mide el tiempo que tarda e imprime el tiempo total y los tiempos de 
        cada par. Usa los mismos parámetros que nuestra implementación para una comparación justa.

        ------------
        Return:
            tiempos_pares: Diccionario que contiene los tiempos de cada par.
            tiempo_total: Tiempo total de ejecución.
        ------------
        """
 
        print("Medicion de tiempo: Smith-Waterman BioPython.")

        # Diccionario donde guardamos el tiempo de cada par.
        tiempos_pares = {}
 
        #Verificamos que hay al menos 2 secuencias
        if len(self.sequences) < 2:
            print("Error: Se necesitan al menos 2 secuencias para comparar.\n")
            return {}, 0
 
        # Obtenemso todas las combinaciones de pares de secuencias sin repetir
        pares_secuencias = list(combinations(range(len(self.sequences)), 2))
        print(f"Número total de pares a alinear: {len(pares_secuencias)}.\n")
 
        #Configuramos primero el alineador de BioPython con los mismos parámetros
        alineador = PairwiseAligner()
        alineador.substitution_matrix = substitution_matrices.load(self.config["matrix_name"])
        alineador.open_gap_score = self.config["gap_penalty"]
        alineador.extend_gap_score = self.config["gap_penalty"]
        alineador.mode = 'local'  #smith-waterman es alineamiento local

        # Registramos tiempo de inicio
        tiempo_inicio_total = time.time()

        #Ejecutamos alineador de Biopython para cada par.
        for indice_seq1, indice_seq2 in pares_secuencias:
            seq1 = str(self.sequences[indice_seq1].seq)
            seq2 = str(self.sequences[indice_seq2].seq)

            # Tiempo inicio del par.
            inicio_par = time.time()

            # Ejecutamos todo el alineamiento completo.
            alineamientos = alineador.align(seq1, seq2)
            alineamientos[0]

             # Tiempo fin del par.
            fin_par = time.time()

            # Tiempo total del par.
            tiempo_par_total = fin_par - inicio_par

            # Actulizamos el diccionario del tiempo de cada par.
            tiempos_pares[(indice_seq1, indice_seq2)] = tiempo_par_total

            print(f"Tiempo del par del algoritmo de biopython ({indice_seq1}, {indice_seq2}) es {tiempo_par_total:.4f} segundos.")

        # Registramos tiempo de fin y calculamos el total
        tiempo_fin_total = time.time()
        tiempo_total = tiempo_fin_total - tiempo_inicio_total
 
        print(f"Tiempo total del algoritmo de biopython: {tiempo_total:.4f} segundos.\n")
 
        return tiempos_pares, tiempo_total
    
    # Octavo paso.
    def comparar_tiempos(self):
        """
        Imprime por pantalla los tiempos de ejecución de ambas versiones de Smith-Waterman y los compara cuál es más rápida.
        """

        print(f"[{self.obtener_hora()}] Paso 8: Comparación de tiempos.\n")

        # Ejecutamos la medición con versión propia:
        # Utilizamos "_" como variable descartada porque la función devuelve dos valores (resultado y tiempo), pero solo nos interesa el tiempo de ejecución.
        _, tiempo_propia = self.medir_tiempo_smith_waterman_propia()

        # Ejecutamos la medición con versión de BioPython:
        _, tiempo_biopython = self.medir_tiempo_smith_waterman_biopython()

        # Imprimimos los resultados:
        print("Resultados:")
        print(f"Smith-Waterman propia: {tiempo_propia:.4f} segundos")
        print(f"Smith-Waterman BioPython: {tiempo_biopython:.4f} segundos.\n")

        # Calculamos cual es la forma más rápida
        if tiempo_biopython < tiempo_propia:
            diferencia = tiempo_propia - tiempo_biopython
            print(f"BioPython es {diferencia:.2f} segundos más rápida.\n")

        elif tiempo_propia < tiempo_biopython:
            diferencia = tiempo_biopython - tiempo_propia
            print(f"Versión propia es {diferencia:.2f} segundos más rápida.\n")
        
        else:
            print("Ambas implementaciones tienen el mismo tiempo de ejecución.\n")

    # Noveno paso.
    def save_sequences(self, output_path, output_format):
        """
        Guarda las secuencias actuales del pipeline en un fichero utilizando SeqIO.write con el formato especificado.

        Tras la escritura, muestra por pantalla un resumen del resultado final, incluyendo el número de secuencias guardadas, el fichero de salida, el formato utilizado 
        y estadísticas básicas como longitud mínima, máxima, media y número de secuencias por tipo biológico.

        ------------
        Params:
            output_path: Ruta de salida.
            output_format: Formato de salida.
        ------------
        """

        print(f"[{self.obtener_hora()}] Paso 9: Escritura de resultados.\n")

        secuencias_escritas = SeqIO.write(self.sequences, output_path, output_format)

        # Resumen del resultado final.
        print(f"Secuencias escritas: {secuencias_escritas}")
        print(f"Fichero de salida: {output_path}")
        print(f"Formato de salida: {output_format}\n")

        print("Almacenamiento de secuencias correcto.\n")
        
    def run (self, ruta):
        """
        Coordina la ejecución completa del pipeline bioinformático.

        Define el flujo de procesamiento llamando de forma ordenada a las distintas etapas: carga de secuencias, procesamiento básico, filtrado, clasificación biológica, 
        cálculo de estadísticas y escritura de resultados.

        Además, muestra mensajes por pantalla que permiten seguir el progreso del pipeline.

        ------------
        Params:
            ruta: Ruta del fichero json.
        ------------
        """

        print("Comienzo del Pipeline:\n")

        self.read_config(ruta) # Lectura del fichero de configuración.
        self.load_sequences(self.config["input_file"], self.config["input_format"]) # Carga de datos.
        self.process() # Procesamiento de los datos.
        self.filter_by_length() # Filtrado por longitud.
        self.classify_and_normalize() # Clasificación y normalización biológica.
        self.compute_basic_stats() # Cálculo de estadísticas finales.

         #Algoritmo local smith-waterman:
        seq1_alineada, seq2_alineada, puntuacion, matchs, mismatchs, gaps = self.smith_waterman(str(self.sequences[0].seq), str(self.sequences[1].seq), self.config["matrix_name"], self.config["gap_penalty"])

        print(f"[{self.obtener_hora()}] Paso 7: Algoritmo local de Smith-Waterman.\n")

        print(f"Score máximo: {puntuacion}")
        print(f"Matchs: {matchs}")
        print(f"Mismatchs: {mismatchs}")
        print(f"Gaps: {gaps}\n")
        print(f"Alineamiento 1: {seq1_alineada}\n")
        print(f"Alineamiento 2: {seq2_alineada}\n")

        print("Algoritmo Smith-Waterman correcto.\n")

        # No imprimimos los tiempos por separado debido a que se recogen en el método comparar_tiempo().
        #self.medir_tiempo_smith_waterman_propia()
        #self.medir_tiempo_smith_waterman_biopython()

        self.comparar_tiempos() # Comparación tiempos del algoritmo Smith-Waterman propio y con Biopython.
        self.save_sequences(self.config["output_file"], self.config["output_format"]) # Guardado de las secuencias.

        print(f"[{self.obtener_hora()}]Finalización del Pipeline.")

"""
Punto de entrada del programa cuando se ejecuta desde línea de comandos.

Este bloque se encarga de:
    1. Validar que se han proporcionado los argumentos necesarios.
    2. Leer los parámetros de entrada (fichero de entrada, ruta de salida y formato de salida).
    3. Crear una instancia de la clase Pipeline.
    4. Ejecutar el pipeline completo mediante el método run().

Permite que el programa funcione como una herramienta ejecutable desde terminal, facilitando su uso en flujos automatizados.
"""

# Ejecución desde línea de comandos.
if __name__ == "__main__":
    if len(sys.argv)<2:
        print("Uso:")
        print("python Pipeline.py SECUENCIA")
        sys.exit(1)

    ruta = sys.argv[1] # Fichero de configuración.

    pipeline = Pipeline()
    pipeline.run(ruta)
