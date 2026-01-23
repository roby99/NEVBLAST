def windw():
    import tkinter as tk
    from Bio.Blast import NCBIWWW
    window = tk.Tk()
    window.title("BLAST")

    # declaring string variable
    # for storing the matrix type,
    # database, sequence, e-value,
    # hit list length, and file name
    matrix_var = tk.StringVar()
    database_var = tk.StringVar()
    sequence_var = tk.StringVar()
    eval_var = tk.StringVar()
    hitList_var = tk.StringVar()
    file_var = tk.StringVar()
    organism_var = tk.StringVar()

    # Dropdown menu options
    dataOptions = [
        "nr",
        "refseq_select",
        "refseq_protein",
        "landmark",
        "swissprot",
        "pataa",
        "env_nr",
        "tsa_nr",
        "pdb"
    ]

    # Dropdown menu options
    matrixOptions = [
        "PAM30",
        "PAM70",
        "PAM250",
        "BLOSUM80",
        "BLOSUM62",
        "BLOSUM45",
        "BLOSUM50",
        "BLOSUM90",
    ]
    gap = {
        "PAM30": "9 1",
        "PAM70": "10 1",
        "PAM250": "14 2",
        "BLOSUM80": "10 1",
        "BLOSUM62": "11 1",
        "BLOSUM45": "15 2",
        "BLOSUM50": "13 2",
        "BLOSUM90": "10 1"
    }
    # datatype of menu text
    dataclicked = tk.StringVar()
    matrixclicked = tk.StringVar()

    # defining a function that will
    # get the variables and runs the
    # blast with them

    def submit():
        matrix = matrix_var.get()
        database = database_var.get()
        sequence = sequence_var.get()
        eval = eval_var.get()
        hits = hitList_var.get()
        file = file_var.get()
        organism = '"' + organism_var.get() + '" [Organism]'

        print("Blasting " + sequence + " using the " + matrixclicked.get() + " matrix and "
              + dataclicked.get() + " database. The e-values are no more than " + eval +
              " and there are " + hits + " hits. The file will be stored in " +
              file + ".xml.")
        window.destroy()
        if organism == '"" [Organism]':
            result_handle = NCBIWWW.qblast("blastp", dataclicked.get(), sequence=sequence, expect=float(eval),
                                           hitlist_size=int(hits), matrix_name=matrixclicked.get(), gapcosts=gap.get(matrixclicked.get()))
        else:
            result_handle = NCBIWWW.qblast("blastp", database, sequence=sequence, expect=float(eval),
                                           hitlist_size=int(hits), matrix_name=matrix, entrez_query=organism)
        with open((file + ".xml"), "w") as out_handle:
            out_handle.write(result_handle.read())
        result_handle.close()

        print("done")

    # created the matrix entry line
    matrix_label = tk.Label(window, text='Matrix', font=('calibre', 10, 'bold'))
    matrix_entry = tk.OptionMenu(window, matrixclicked, *matrixOptions)

    # created the database entry line
    database_label = tk.Label(window, text='Database', font=('calibre', 10, 'bold'))
    database_entry = tk.OptionMenu(window, dataclicked, *dataOptions)

    # created the sequence entry line
    sequence_label = tk.Label(window, text='Sequence or Accessions Number', font=('calibre', 10, 'bold'))
    sequence_entry = tk.Entry(window, textvariable=sequence_var, font=('calibre', 10, 'normal'))

    # created the e-value entry line
    eval_label = tk.Label(window, text='E-value', font=('calibre', 10, 'bold'))
    eval_entry = tk.Entry(window, textvariable=eval_var, font=('calibre', 10, 'normal'))

    # created the hits entry line
    hits_label = tk.Label(window, text='Number of Hits', font=('calibre', 10, 'bold'))
    hits_entry = tk.Entry(window, textvariable=hitList_var, font=('calibre', 10, 'normal'))

    # created the file name entry line
    file_label = tk.Label(window, text='File Name', font=('calibre', 10, 'bold'))
    file_entry = tk.Entry(window, textvariable=file_var, font=('calibre', 10, 'normal'))

    # created the file name entry line
    organism_label = tk.Label(window, text='Organism', font=('calibre', 10, 'bold'))
    organism_entry = tk.Entry(window, textvariable=organism_var, font=('calibre', 10, 'normal'))

    # creating a button using the widget
    # Button that will call the submit function
    sub_btn = tk.Button(window, text='Submit', command=submit)

    # placing the label and entry in
    # the required position using grid
    # method
    matrix_label.grid(row=0, column=0)
    matrix_entry.grid(row=0, column=1)
    database_label.grid(row=1, column=0)
    database_entry.grid(row=1, column=1)
    sequence_label.grid(row=2, column=0)
    sequence_entry.grid(row=2, column=1)
    eval_label.grid(row=3, column=0)
    eval_entry.grid(row=3, column=1)
    hits_label.grid(row=4, column=0)
    hits_entry.grid(row=4, column=1)
    file_label.grid(row=5, column=0)
    file_entry.grid(row=5, column=1)
    organism_label.grid(row=6, column=0)
    organism_entry.grid(row=6, column=1)
    sub_btn.grid(row=7, column=1)
    window.mainloop()

