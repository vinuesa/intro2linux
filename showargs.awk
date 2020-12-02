# showargs.awk v0.3
# bloque BEGIN{} de inicializacion
BEGIN {
   RS=">"
   FS="\n"
   printf "BEGIN block: A=%d, B=%s, NF=%s, NR=%s, FNR=%s\n", A, B, NF, NR, FNR
   for (i=0; i < ARGC; i++)
      printf "\tARGV[%d] = %s\n", i, ARGV[i] 
   print "End of BEGIN block\n--------------------------------\n"
}

# filtro
NR > 1 && $1 == B

# bloque END{} de procesamiento final
END { printf "END block: A=%d, B=%s, NF=%s, NR=%s, FNR=%s\n", A, B, NF, NR, FNR }
END { for (i=1; i<=NF; i++) print $i }
