 # inicialización para leer archivos fasta
BEGIN{RS=">";FS="\n"}

# cuando trabajamos con el primer archivo de la lista de argumentos, 
# que contiene la lista de cabeceras FASTA a filtrar (formato: >etiqueta1\n>etiqueta2 ...)
#  se cumple que NR==FNR y simplemente inicializamos el hash labels_h
NR == FNR { labels_h[$1]++ }

# NR>FNR se cumple cuando awk comienza a procesar el segundo archivo,
#  el FASTA con las secuencias a filtrar
#  Si la etiqueta de la cabecera FASTA existe en el hash_
NR > FNR { if ($1 in labels_h && $0!="") printf ">%s", $0 }
