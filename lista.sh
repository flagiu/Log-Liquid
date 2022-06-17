# elenca le configurazioni in ordine di indice, tuttavia
# ls include anche file ../Cnf*-esteso.tpr ed eventuali file zippati.

# il comando awk esclude il file .tpr
# lo script in python esclude tutti i files che non rispettano il formato
# ../Cnf*-{intero} (quindi si potrebbe eliminare awk)
# inoltre produce lista.dat (nome,indice,tempo) e listatempi.dat (indice,tempo)

ls ../Cnf*-* | sort -n -t "-" -k 2 | awk '{if (NR>1) print $0}' > listaconf.dat
python3 lista.py
rm listaconf.dat
