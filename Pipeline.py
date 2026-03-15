import sys
from datetime import datetime
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import UndefinedSequenceError


# Ejercicio 2: Estado interno del pipeline.
# Mantiene un estado interno persistente que permite alamacenar datos y resultados intermedios.
class Pipeline:
    """
    Clase que define la arquitectura básica de un pipeline bioinformático.
    """

    def __init__(self, input_path, input_format):
        """
        Inicialización estado interno del pipeline. 

        ---------
        Params:
        input_path: str
           Ruta al fichero donde se encuentran las secuencias.
        input_format: str
            Formato del fichero. 
        ---------
        """

        self.input_path = input_path
        self.input_format = input_format # Guardamos nombre y formato del fichero

        # Diccionario que almacenará las secuencias como objetos.
        self.sequences = {}

        # Diccionario que guarda información descriptiva sobre secuencias.
        self.metadata = {
            "n_dna":0,
            "n_rna":0,
            "n_prot":0,
            "n_total":0,
            "n_min":0,
            "n_max":0,
            "n_mean":0,
            "n_seq":0

        }

        # Diccionario con los parámetros de configuración del pipeline.
        self.config = {
            "min_length" : 5  # Ponemos como parámetro de filtrado por longitud 5
        }

        # Diccionario con los resultados finales.
        self.results = {}

    # Ejercicio 7: Carga de secuencias con SeqIO
    def load_sequences(self):
        """
        Carga secuencias desde un fichero usando SeqIO.parse y las almacena en un diccionario indexado por ID.
        """
        print(f"[{datetime.now().strftime('%H:%M:%S')}]Paso 1: Cargando secuencias...")

        for record in SeqIO.parse(self.input_path, self.input_format):
            self.sequences[record.id] = record  # Para poder acceder por ID fácilemnte.
        
        # Ante fichero vacío o incorrecto: 
        if len(self.sequences) == 0:
            print(f"[{datetime.now().strftime('%H:%M:%S')}]No se cargaron secuencias.")
            return
        
        longitudes =  [len(record.seq) for record in self.sequences.values()]

        self.metadata["n_total"] = len(self.sequences)
        self.metadata["n_min"] = min(longitudes)
        self.metadata["n_max"] = max(longitudes)
        self.metadata["n_mean"] = sum(longitudes) / len(longitudes)

        print(f"[{datetime.now().strftime('%H:%M:%S')}]Se cargaron {self.metadata['n_total']} secuencias.")

    # Ejercicio 8: Filtrado por longitud:
    def filter_by_length(self, min_length):
        """
        Elimina secuencias de longitud inferior a un umbral establecido anteriormente (5).
        """
        print(f"[{datetime.now().strftime('%H:%M:%S')}]Paso 2: Filtrando por longitud mínima...")

        filtradas = {
            id_: record
            for id_, record in self.sequences.items()
            if len(record.seq) >= min_length
        }

        self.sequences = filtradas

        print(f"[{datetime.now().strftime('%H:%M:%S')}]Secuencias resultantes tras el filtrado: {len(self.sequences)}")

        # Actualizamos el metadata tras filtrado:
        if len(self.sequences) > 0:
            longitudes = [len(record.seq) for record in self.sequences.values()]
            self.metadata["n_total"] = len(self.sequences)
            self.metadata["n_min"] = min(longitudes)
            self.metadata["n_max"] = max(longitudes)
            self.metadata["n_mean"] = sum(longitudes) / len(longitudes)

    #Ejercicio 5: Calsificación y normalización biológica de secuencias.
    def classify_and_normalize(self) :
        """
        Clasifica las secuencias y convierte ARN --> ADN. Descarta secuencias inválidas.
        """
        print(f"[{datetime.now().strftime('%H:%M:%S')}]Paso 3: Clasificación y normalización secuencias...")

        dict_dna = {"A", "C", "G", "T"}
        dict_rna = {"A", "C", "G", "U"}
        dict_prot = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"}

        seq_valida = {}

        for record_id, record in self.sequences.items():
            try:
                seq = str(record.seq)
            except UndefinedSequenceError:
                print(f"[{datetime.now().strftime('%H:%M:%S')}]Secuencia sin contenido: {record_id}")
                continue

            letras = set(seq)
            
            # Comprobación es ADN
            if letras.issubset(dict_dna):
                self.metadata["n_dna"] += 1
                seq_valida[record_id] = record

            # Comprobación es ARN
            elif letras.issubset(dict_rna):
                self.metadata["n_rna"] += 1
                dna_str = str(record.seq).replace("U", "T") # Reemplaza U por T
                record.seq = Seq(dna_str)
                seq_valida[record_id] = record
            
            # Comprobación es proteína
            elif letras.issubset(dict_prot):
                self.metadata["n_prot"] += 1
                seq_valida[record_id] = record
            
            # Si no es ninguna comprobación de las anteriores
            else:
                print(f"[{datetime.now().strftime('%H:%M:%S')}]Secuencia inválida: {record_id}")

        self.sequences = seq_valida
        print(f"[{datetime.now().strftime('%H:%M:%S')}]Clasificación completada.")

    #Ejercicio 4: Soporte para múltiples secuencias por línea de comandos.
    def process(self):
        """
        Realiza un procesamiento simple sobre la secuencia
        """
        print(f"[{datetime.now().strftime('%H:%M:%S')}]Paso 4: Procesando secuencias...")

        longitud = [len(record.seq) for record in self.sequences.values()] # Calcular longitud
        self.results["longitud"] = longitud
        print(f"[{datetime.now().strftime('%H:%M:%S')}]Longitud de cada cadena: {longitud}")

    #Ejercicio 6: Estadísticas globales del pipeline.
    def compute_basic_stats(self):
        """
        Calcula estadísticas descriptivas globales sobre las secuencias actualmente almacenadas.
        """
        print(f"[{datetime.now().strftime('%H:%M:%S')}]Paso 5: Calculando estadísticas globales...")

        total = len(self.sequences)
        if total == 0:
            print(f"[{datetime.now().strftime('%H:%M:%S')}]No hay secuencias válidas.")
            return
        
        longitud = [len(record.seq) for record in self.sequences.values()]

        min_long = min(longitud)
        max_long = max(longitud)
        media_long = sum(longitud) / total

        self.metadata["total_secuencias"] = total
        self.metadata["min_longitud"] = min_long
        self.metadata["max_longitud"] = max_long
        self.metadata["media_longitud"] = media_long

        # Número total de secuencias.
        print(f"[{datetime.now().strftime('%H:%M:%S')}]Número total de secuencias: {total}") 

        # Longitud mínima, máxima y media.
        print(f"[{datetime.now().strftime('%H:%M:%S')}]Longitud mínima: {min_long}")
        print(f"[{datetime.now().strftime('%H:%M:%S')}]Longitud máxima: {max_long}")
        print(f"[{datetime.now().strftime('%H:%M:%S')}]Longitud media: {media_long}")
        print(f"[{datetime.now().strftime('%H:%M:%S')}]Número de secuencias de ADN: {self.metadata['n_dna']}")
        print(f"[{datetime.now().strftime('%H:%M:%S')}]Número de secuencias de ARN: {self.metadata['n_rna']}")
        print(f"[{datetime.now().strftime('%H:%M:%S')}]Número de secuencias de proteína: {self.metadata['n_prot']}")

        # Número de secuencias de cada tipo biológico.

    # Ejercicio 9: Escritura de resultados:
    def save_sequences(self, output_path, output_format):
        """
        Guarda las secuencias actuales usando SeqIO.write.
        """
        print(f"[{datetime.now().strftime('%H:%M:%S')}]Paso 6: Guardando secuencias...")

        SeqIO.write(self.sequences.values(), output_path, output_format)

        print(f"[{datetime.now().strftime('%H:%M:%S')}]Fichero generado: {output_path}")
        print(f"[{datetime.now().strftime('%H:%M:%S')}]Número final de secuencias: {len(self.sequences)}")

    
    def run(self):
        """
        Punto de entrada principal del pipeline. Controla el flujo de ejecución.
        """
        print(f"[{datetime.now().strftime('%H:%M:%S')}]Inicio del Pipeline")

        self.load_sequences()
        self.filter_by_length(self.config["min_length"])
        self.classify_and_normalize()
        self.process()
        self.compute_basic_stats()
        self.save_sequences("output.fasta", self.input_format)

        print(f"[{datetime.now().strftime('%H:%M:%S')}]Fin del Pipeline")
    
   



# Ejecución desde línea de comandos.
if __name__ == "__main__":
    if len(sys.argv)<3:
        print("Uso:")
        print("python Pipeline.py SEQ1")
        sys.exit(1)

    input_path = sys.argv[1]     # fichero de secuencias
    input_format = sys.argv[2]   # formato del fichero

    pipeline = Pipeline(input_path, input_format)
    pipeline.run()


        
