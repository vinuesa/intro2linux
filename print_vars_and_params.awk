# print_vars_and_params.awk imprime variables y parametros pasados al script como se muestra:
# USO:  awk -v var1=$RANDOM -f ./print_vars_and_params.awk param1=valParam1 file1
BEGIN { 
    print "# ARGV contiene los siguientes argumentos posicionales:"
    
    # recorremos los elementos de ARGV, usando como indices posicionales los guardados en ARGC
    #  y los imprimimos
    for (i=0; i < ARGC; i++) print ARGV[i]
    
    # imprimimos el valor de ARGC
    print "\n# ARGC == lenght(ARGV) =>", length(ARGV), "\nARGC =", ARGC
    
    # imprimimos los valores de los parametros y variables pasados al script
    #   noten el uso de "[" y "]" para facilitar la visualizacion del contenido de cada param o var
    print "["param1"]", "["param2"]", "["var1"]", "["var2"]"
}
