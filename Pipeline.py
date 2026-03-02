import sys
from Bio.Seq import Seq

# Ejercicio 2: Estado interno del pipeline.
# Mantiene un estado interno persistente que permite alamacenar datos y resultados intermedios.
class Pipeline:
    """
    Clase que define la arquitectura básica de un pipline bioinformático.
    """

    def __init__(self, input_sequence):
        """
        Inicialización estado interno del pipeline. 

        ---------
        Params:
        input_sequence: list
           Lista de secuencias pasadas por línea de comandos. 
        ---------
        """

        self.input_sequence = input_sequence # Secuencia cruda.

        # Lista que almacenará las secuencias como objetos.
        self.sequences = []

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

        # Diccionario con los parámetros de configuración del pipiline.
        self.config = {}

        # Diccionario con los resultados finales.
        self.results = {}

    #Ejercicio 3: Flujo de ejecución. Etapas del pipeline.
    # Ejercicio 7: Carga de secuencias con SeqIO
    def load_sequences(self):
        """
        Carga secuencias desde un fichero usando SeqIO.parse y las almacena en un diccionario indexado por ID.
        """
        print("Paso 1: Cargando secuencias...")

        for seq in SeqIO.parse():
            self.sequences.append(Seq(seq.upper()))

        print(f"Se cargaron {len(self.sequences)} secuencias")

    #Ejercicio 5: Calsificación y normalización biológica de secuencias.
    def classify_and_normalize(self) :
        """
        Calsifica las secuencias y convierte ARN --> ADN. Descarta secuencias inválidas.
        """
        print("Paso 2: Clasificación y normalización secuencias...")

        dict_dna = {"A", "C", "G", "T"}
        dict_rna = {"A", "C", "G", "U"}
        dict_prot = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"}

        seq_valida = []

        for seq in self.sequences:
            letras = set(str(seq))

            # Comprobación es ADN
            if letras.issubset(dict_dna):
                self.metadata["n_dna"] += 1
                seq_valida.append(seq)

            # Comprobación es ARN
            elif letras.issubset(dict_rna):
                self.metadata["n_rna"] += 1
                dna_str = str(seq).replace("U", "T") # Reemplaza U por T
                dna_seq = Seq(dna_str)
                seq_valida.append(dna_seq)
            
            # Comprobación es proteína
            elif letras.issubset(dict_prot):
                self.metadata["n_prot"] += 1
                seq_valida.append(seq)
            
            # Si no es ninguna comprobación de las anteriores
            else:
                print(f"Secuencia inválida: {seq}")

        self.sequences = seq_valida
        print("Clasificación completada.")

    #Ejercicio 4: Soporte para múltiples secuencias por línea de comandos.
    def process(self):
        """
        Realiza un procesamiento simple sobre la secuencia
        """
        print("Paso 3: Procesando secuencias...")

        longitud = [len(seq) for seq in self.sequences] # Calcular longitud
        self.results["longitud"] = longitud
        print(f"Longitud de cada cadena: {longitud}")

    #Ejercicio 6: Estadísticas globales del pipeline.
    def compute_basic_stats(self):
        """
        Calcula estadísticas descriptivas globales sobre las secuencias actualmente almacenadas.
        """
        print("Paso 4: Calculando estadísticas globales...")

        total = len(self.sequences)
        if total == 0:
            print("No hay secuencias válidas")
            return
        
        longitud = [len(seq) for seq in self.sequences]

        min_long = min(longitud)
        max_long = max(longitud)
        media_long = sum(longitud) / total

        self.metadata["total_secuencias"] = total
        self.metadata["min_longitud"] = min_long
        self.metadata["max_longitud"] = max_long
        self.metadata["media_longitud"] = media_long

        # Número total de secuencias.
        print(f"Número total de secuencias:{total}") 

        # Longitud mínima, máxima y media.
        print(f"Longitud mínima: {min_long}")
        print(f"Longitud máxima: {max_long}")
        print(f"Longitud media: {media_long}")
        print(f"Número de secuencias de ADN: {self.metadata['n_dna']}")
        print(f"Número de secuencias de ARN: {self.metadata['n_rna']}")
        print(f"Número de secuencias de proteína: {self.metadata['n_prot']}")

        # Número de secuencias de cada tipo biológico.
    
    def run(self):
        """
        Punto de entrada principal del pipeline. Controla el flujo de ejecución.
        """
        print("Inicio del Pipeline")

        self.load()
        self.classify_and_normalize()
        self.process()
        self.compute_basic_stats()

        print("Fin del Pipeline")
    
   



# Ejecución desde línea de comandos.
if __name__ == "__main__":
    if len(sys.argv)<2:
        print("Uso:")
        print("python Pipline.py SEQ1")
        sys.exit(1)

    input_sequence = sys.argv[1:] # Secuencias como listas

    pipeline = Pipeline(input_sequence)
    pipeline.run()


        
