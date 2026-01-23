import tkinter as tk
from tkinter import ttk, messagebox
import threading
import time
from PIL import Image, ImageTk
import NEVBlastBlast as nbb
import NEVBlast_Backend as nev
import Plot
from Bio import Entrez
Entrez.email = "nev-blast@uwp.edu"  # Required by NCBI

window = tk.Tk()
window.title("NEV-BLAST: Novel Variant BLAST Analysis Tool")
window.geometry('1400x850')
window.configure(bg='#f5f5f5')
print("*** RUNNING THE CORRECT FILE - VERSION 2024-11-08 ***")

file_var = tk.StringVar()
sig_var = tk.StringVar()
sequence_var = tk.StringVar(value="AAM07687.1")
matrix_var = tk.StringVar()
database_var = tk.StringVar()
eval_var = tk.StringVar(value="0.0001")
hits_var = tk.StringVar(value="100")
organism_var = tk.StringVar()
cancel_requested = False

# UW-Parkside colors
UWP_GREEN = '#1e5631'
UWP_LIGHT_GREEN = '#2d7a4a'
UWP_WHITE = '#ffffff'
UWP_BLACK = '#2d2d2d'
UWP_GRAY = '#f5f5f5'

# Create main container
main_container = tk.Frame(window, bg=UWP_GRAY)
main_container.pack(fill='both', expand=True, padx=25, pady=(25, 0))

# LEFT PANEL
left_panel = tk.Frame(main_container, bg=UWP_WHITE, relief='solid', borderwidth=2)
left_panel.grid(row=0, column=0, sticky='nsew', padx=(0, 15))

# RIGHT PANEL
right_panel = tk.Frame(main_container, bg=UWP_WHITE, relief='solid', borderwidth=2)
right_panel.grid(row=0, column=1, sticky='nsew')

# Configure grid weights
main_container.columnconfigure(0, weight=1, minsize=650)
main_container.columnconfigure(1, weight=1, minsize=650)
main_container.rowconfigure(0, weight=1)

# ============= LEFT PANEL CONTENT =============

# Header section with logo and title side by side
header_frame = tk.Frame(left_panel, bg=UWP_WHITE)
header_frame.pack(pady=(10, 5), padx=(0))

# Logo on the left
logo_frame = tk.Frame(header_frame, bg=UWP_WHITE)
logo_frame.pack(side='left', padx=(0, 0), anchor='n')

try:
    logo_image = Image.open("uwp_logo.png")
    logo_image = logo_image.resize((300, 100), Image.Resampling.LANCZOS)
    logo_photo = ImageTk.PhotoImage(logo_image)
    logo_label = tk.Label(logo_frame, image=logo_photo, bg=UWP_WHITE)
    logo_label.image = logo_photo
    logo_label.pack(padx=0, pady=0)

except:
    tk.Label(logo_frame, text="UW-PARKSIDE", 
            font=('Arial', 18, 'bold'), bg=UWP_WHITE, fg=UWP_GREEN).pack()

# Title on the right of logo
title_frame = tk.Frame(header_frame, bg=UWP_WHITE)
title_frame.pack(side='left', anchor='n', pady=(3, 0))

tk.Label(title_frame, text="NEV-BLAST", 
         font=('Arial', 40, 'bold'), bg=UWP_WHITE, fg=UWP_GREEN).pack(anchor='w')
tk.Label(title_frame, text="Novel Variant BLAST Analysis", 
         font=('Arial', 14, 'italic'), bg=UWP_WHITE, fg=UWP_BLACK).pack(anchor='w')

# Separator line
tk.Frame(left_panel, height=2, bg=UWP_GREEN).pack(fill='x', padx=30, pady=5)


# What is NEV-BLAST section - cleaner design
desc_container = tk.Frame(left_panel, bg=UWP_WHITE)
desc_container.pack(fill='x', padx=20, pady=8)

tk.Label(desc_container, text="What is NEV-BLAST?", 
         font=('Arial', 20, 'bold'), bg=UWP_WHITE, fg=UWP_GREEN).pack(anchor='w', pady=(0, 6))

desc_frame = tk.Frame(desc_container, bg=UWP_WHITE, relief='solid', borderwidth=2)
desc_frame.pack(fill='x')

desc_inner = tk.Frame(desc_frame, bg=UWP_WHITE)
desc_inner.pack(fill='both', padx=15, pady=12)

tk.Label(desc_inner, 
        text="NEV-BLAST identifies protein variants by:\n"
             "• Searching protein databases for similar sequences\n"
             "• Analyzing signature amino acid positions\n"
             "• Scoring conservation using BLOSUM matrices\n"
             "• Visualizing variants in 3D interactive plots",
        font=('Arial', 16), bg=UWP_WHITE, justify='left', wraplength=550).pack(anchor='w')


# How to Use section - cleaner box design
instr_container = tk.Frame(left_panel, bg=UWP_WHITE)
instr_container.pack(fill='x', padx=20, pady=8)

# Title bar
tk.Label(instr_container, text="How to Use", 
         font=('Arial', 20, 'bold'), bg=UWP_WHITE, fg=UWP_GREEN).pack(anchor='w', pady=(0, 6))

# Content box with border
instr_frame = tk.Frame(instr_container, bg=UWP_WHITE, relief='solid', borderwidth=2)
instr_frame.pack(fill='x')

instr_inner = tk.Frame(instr_frame, bg=UWP_WHITE)
instr_inner.pack(fill='both', padx=15, pady=12)

tk.Label(instr_inner,
        text="1. Enter a unique file name for your analysis\n"
             "2. Define signature sequences as: [['position','amino acid'], ...]\n"
             "    - An example would be [['5','C'],['54','C']]\n"
             "3. Enter protein sequence or NCBI accession number\n"
             "4. Adjust BLAST parameters (optional)\n"
             "5. Click 'RUN COMPLETE ANALYSIS' button\n"
             "6. Wait approximately 15-20 minutes for completion",
        font=('Arial', 16), bg=UWP_WHITE, justify='left', wraplength=550).pack(anchor='w')

# Expected Output section - cleaner box design
output_container = tk.Frame(left_panel, bg=UWP_WHITE)
output_container.pack(fill='x', padx=20, pady=8)

# Title bar
tk.Label(output_container, text="Expected Output", 
         font=('Arial', 20, 'bold'), bg=UWP_WHITE, fg=UWP_GREEN).pack(anchor='w', pady=(0, 6))

# Content box with border
output_frame = tk.Frame(output_container, bg=UWP_WHITE, relief='solid', borderwidth=2)
output_frame.pack(fill='x')

output_inner = tk.Frame(output_frame, bg=UWP_WHITE)
output_inner.pack(fill='both', padx=15, pady=12)

tk.Label(output_inner,
        text="1. filename.xml - Raw BLAST search results from NCBI\n"
             "2. filename.txt - Parsed signature data with position mapping\n"
             "3. filename.csv - Conservation scores for each variant\n"
             "4. Interactive 3D Plot - Browser visualization of signature conservation across\n"
             "    E-values and variants",
        font=('Arial', 16), bg=UWP_WHITE, justify='left', wraplength=550).pack(anchor='w')

# Version info
tk.Label(left_panel, text="Version 3.0 | Optimized for 1000+ sequences | © 2025 UW-Parkside", 
        font=('Arial', 12, 'bold'), bg=UWP_WHITE, fg="#070707").pack(side='bottom', pady=8)


# ============= RIGHT PANEL CONTENT =============

# REQUIRED INPUTS - Green bar at top of right panel (outside scrollable area)
required_header = tk.Frame(right_panel, bg=UWP_GREEN, bd=0)
required_header.pack(fill='x', side='top')
tk.Label(required_header, text="REQUIRED INPUTS", 
         font=('Arial', 26, 'bold'), bg=UWP_GREEN, fg=UWP_WHITE).pack(pady=8)

# Scrollable canvas for content
canvas = tk.Canvas(right_panel, bg=UWP_WHITE, highlightthickness=0)
scrollable_frame = tk.Frame(canvas, bg=UWP_WHITE)

scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
canvas_window = canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")

def on_canvas_configure(event):
    canvas.itemconfig(canvas_window, width=event.width)
canvas.bind("<Configure>", on_canvas_configure)


# ============= ORGANISM VALIDATION FUNCTIONS =============

organism_cache = {}  # Cache search results
      
def search_organism(query):
    """Search NCBI Taxonomy for organism names with taxid - shows multiple results"""
    print(f"[DEBUG] Searching for: '{query}'")
    
    if not query or len(query) < 3:
        print("[DEBUG] Query too short")
        return []
    
    # Check cache first
    if query.lower() in organism_cache:
        print(f"[DEBUG] Found in cache: {len(organism_cache[query.lower()])} results")
        return organism_cache[query.lower()]
    
    try:
        # FIXED: Add [Organism] field tag to prevent 400 Bad Request error
        search_term = f"{query}[Organism]"
        print(f"[DEBUG] Calling Entrez.esearch with term: '{search_term}'...")
        
        handle = Entrez.esearch(db="taxonomy", term=search_term, retmax=20)
        record = Entrez.read(handle)
        handle.close()
        
        print(f"[DEBUG] Found {len(record.get('IdList', []))} IDs")
        
        if record["IdList"]:
            # Fetch details for found IDs
            id_list = record["IdList"]
            print(f"[DEBUG] Fetching details for {len(id_list)} organisms...")
            handle = Entrez.efetch(db="taxonomy", id=",".join(id_list), retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            
            # Extract organism names with taxid
            organisms = []
            for rec in records:
                if "ScientificName" in rec and "TaxId" in rec:
                    taxid = rec["TaxId"]
                    name = rec["ScientificName"]
                    # Format: "Organism name (taxid:12345)"
                    display = f"{name} (taxid:{taxid})"
                    organisms.append(display)
                    print(f"[DEBUG] Added: {display}")
            
            # Sort alphabetically for easier selection
            organisms.sort()
            
            organism_cache[query.lower()] = organisms
            print(f"[DEBUG] Returning {len(organisms)} results")
            return organisms
    except Exception as e:
        print(f"[ERROR] Organism search error: {e}")
        import traceback
        traceback.print_exc()
        return []
    
    return []


# Debounce timer to avoid rate limiting
search_timer = None

def update_organism_dropdown(event):
    """Update dropdown with organism suggestions as user types - with debouncing"""
    global search_timer
    
    query = organism_var.get()
    print(f"[DEBUG] KeyRelease event - current text: '{query}'")
    
    # Cancel any pending search
    if search_timer is not None:
        window.after_cancel(search_timer)
        search_timer = None
        print("[DEBUG] Cancelled pending search")
    
    # Clear dropdown if query too short
    if len(query) < 3:
        organism_dropdown['values'] = []
        print("[DEBUG] Query too short, clearing dropdown")
        return
    
    # Check cache first - instant results!
    if query.lower() in organism_cache:
        cached_results = organism_cache[query.lower()]
        organism_dropdown['values'] = cached_results
        print(f"[DEBUG] Using cached results: {len(cached_results)} items")
        
        # Auto-open the dropdown to show cached results
        try:
            organism_dropdown.tk.call('ttk::combobox::Post', organism_dropdown)
            print("[DEBUG] Dropdown opened (from cache)")
        except:
            pass
        return
    
    # Show "Searching..." message immediately
    organism_dropdown['values'] = ["Searching..."]
    print("[DEBUG] Waiting for user to stop typing...")
    
    # Wait 2 seconds after user stops typing before searching
    def delayed_search():
        global search_timer
        search_timer = None
        print(f"[DEBUG] User stopped typing, searching for: '{query}'")
        
        organism_dropdown['values'] = ["Searching NCBI..."]
        
        # Search for organisms in background thread
        def do_search():
            print("[DEBUG] Thread started")
            results = search_organism(query)
            print(f"[DEBUG] Thread finished, got {len(results)} results")
            
            # Update dropdown on main thread
            def update_ui():
                if results:
                    organism_dropdown['values'] = results
                    print(f"[DEBUG] Updated dropdown with {len(results)} results")
                    
                    # Auto-open the dropdown to show results
                    try:
                        organism_dropdown.tk.call('ttk::combobox::Post', organism_dropdown)
                        print("[DEBUG] Dropdown opened automatically")
                    except Exception as e:
                        print(f"[DEBUG] Could not auto-open dropdown: {e}")
                else:
                    organism_dropdown['values'] = ["No results found - try full genus/species name"]
                    print("[DEBUG] No results found")
            
            window.after(0, update_ui)
        
        thread = threading.Thread(target=do_search, daemon=True)
        thread.start()
    
    # Schedule search for 2 seconds from now
    search_timer = window.after(2000, delayed_search)
    print("[DEBUG] Search scheduled in 2 seconds")


# Required inputs content
required_content = tk.Frame(scrollable_frame, bg=UWP_WHITE)
required_content.pack(fill='x', padx=60, pady=8)

# File Name
tk.Label(required_content, text='File Name:', font=('Arial', 18, 'bold'), 
         bg=UWP_WHITE, fg=UWP_BLACK).pack(pady=(0, 4))
file_entry = tk.Entry(required_content, textvariable=file_var, 
                     font=('Arial', 18), width=50, relief='solid', borderwidth=2,
                     justify='center')
file_entry.pack(pady=(0, 4))

# Signature Sequences
tk.Label(required_content, text='Signature Sequences:', font=('Arial', 18, 'bold'), 
         bg=UWP_WHITE, fg=UWP_BLACK).pack(pady=(0, 4))
tk.Label(required_content, text='Format: [[\'5\',\'C\'],[\'54\',\'C\']],[[\'67\',\'S\']]', 
         font=('Arial', 14, 'italic'), bg=UWP_WHITE, fg="#101010").pack(pady=(0, 4))
sig_entry = tk.Entry(required_content, textvariable=sig_var, 
                    font=('Arial', 18), width=50, relief='solid', borderwidth=2,
                    justify='center')
sig_entry.pack(pady=(0, 4))

# Sequence
tk.Label(required_content, text='Protein Sequence or Accession Number:', 
         font=('Arial', 18, 'bold'), bg=UWP_WHITE, fg=UWP_BLACK).pack(pady=(0, 4))
sequence_entry = tk.Entry(required_content, textvariable=sequence_var, 
                         font=('Arial', 18), width=50, relief='solid', borderwidth=2,
                         justify='center')
sequence_entry.pack(pady=(0, 4))

# BLAST PARAMETERS - Green bar inside scrollable area
blast_header = tk.Frame(scrollable_frame, bg=UWP_GREEN, bd=0, highlightthickness=0)
blast_header.pack(fill='x', padx=0, pady=(10, 0), expand=False)
tk.Label(blast_header, text="BLAST PARAMETERS", 
         font=('Arial', 26, 'bold'), bg=UWP_GREEN, fg=UWP_WHITE).pack(pady=8)

# BLAST content
blast_content = tk.Frame(scrollable_frame, bg=UWP_WHITE)
blast_content.pack(fill='x', padx=60, pady=8)

# Matrix
tk.Label(blast_content, text='Scoring Matrix:', font=('Arial', 18, 'bold'), 
         bg=UWP_WHITE, fg=UWP_BLACK).pack(pady=(0, 4))
matrix_options = ["BLOSUM62", "BLOSUM45", "BLOSUM50", "BLOSUM80", "BLOSUM90", "PAM30", "PAM70", "PAM250"]
matrix_dropdown = ttk.Combobox(blast_content, textvariable=matrix_var, values=matrix_options,
                              font=('Arial', 18), state='readonly', width=48, justify='center')
matrix_dropdown.set("BLOSUM62")
matrix_dropdown.pack(pady=(0, 4))

# Database
tk.Label(blast_content, text='Database:', font=('Arial', 18, 'bold'), 
         bg=UWP_WHITE, fg=UWP_BLACK).pack(pady=(0, 4))
database_options = ["refseq_protein", "refseq_select", "nr", "swissprot", "pdb"]
database_dropdown = ttk.Combobox(blast_content, textvariable=database_var, values=database_options,
                                font=('Arial', 18), state='readonly', width=48, justify='center')
database_dropdown.set("refseq_protein")
database_dropdown.pack(pady=(0, 4))

# E-value
tk.Label(blast_content, text='E-value Threshold:', font=('Arial', 18, 'bold'), 
         bg=UWP_WHITE, fg=UWP_BLACK).pack(pady=(0, 4))
eval_entry = tk.Entry(blast_content, textvariable=eval_var, 
                     font=('Arial', 18), width=50, relief='solid', borderwidth=2,
                     justify='center')
eval_entry.pack(pady=(0,4))

# Number of Hits
tk.Label(blast_content, text='Number of Hits:', font=('Arial', 18, 'bold'), 
         bg=UWP_WHITE, fg=UWP_BLACK).pack(pady=(0, 4))
hits_entry = tk.Entry(blast_content, textvariable=hits_var, 
                     font=('Arial', 18), width=50, relief='solid', borderwidth=2,
                     justify='center')
hits_entry.pack(pady=(0, 4))

# Organism with autocomplete
tk.Label(blast_content, text='Organism (Optional):', font=('Arial', 18, 'bold'), 
         bg=UWP_WHITE, fg=UWP_BLACK).pack(pady=(0, 4))
tk.Label(blast_content, text='Type 4+ characters to search NCBI Taxonomy', 
         font=('Arial', 16, 'bold'), bg=UWP_WHITE, fg="#1e5631").pack(pady=(0, 4))
organism_dropdown = ttk.Combobox(blast_content, textvariable=organism_var,
                                font=('Arial', 18), width=48, justify='center')
organism_dropdown.pack(pady=(0, 10))
organism_dropdown.bind('<KeyRelease>', update_organism_dropdown)

canvas.pack(fill="both", expand=True)

# ============= BOTTOM BUTTON PANEL =============

# Create button container at bottom - flush with window edge
button_container = tk.Frame(window, bg='#f5f5f5')
button_container.pack(side='bottom', fill='x', padx=25, before=main_container)

button_inner = tk.Frame(button_container, bg='#f5f5f5', height=55)
button_inner.pack(fill='x', pady=(10, 10))

# RUN button - Using Label styled as button for full color control
run_all_btn = tk.Label(button_inner, 
                       text='🚀 RUN COMPLETE ANALYSIS',
                       bg='#1e5631',
                       fg='white',
                       font=('Arial', 18, 'bold'),
                       relief='raised',
                       borderwidth=4,
                       cursor='arrow')

run_all_btn.pack(side='left', fill='both', expand=True, padx=0, pady=0)

# CANCEL button - Using Label styled as button
cancel_btn = tk.Label(button_inner,
                      text='✖ CANCEL',
                      bg='#d91E18',
                      fg='white',
                      font=('Arial', 16, 'bold'),
                      relief='raised',
                      borderwidth=4,
                      cursor='hand2')

cancel_btn.pack(side='left', fill='y', padx=(5, 0), pady=0, ipadx=30, ipady=15)

# Initially disable cancel button (light red to show it's cancel)
cancel_btn.config(bg='#EE7674', fg='white', cursor='')

# Click bindings (Labels don't have command parameter, use bind instead)
run_all_btn.bind('<Button-1>', lambda e: run_all() if run_all_btn['cursor'] == 'arrow' else None)
cancel_btn.bind('<Button-1>', lambda e: cancel_workflow() if cancel_btn['cursor'] == 'hand2' else None)

# Hover effects for Labels
def on_enter_run(e):
    if run_all_btn['cursor'] == 'arrow':
        run_all_btn['bg'] = '#16412a'

def on_leave_run(e):
    if run_all_btn['cursor'] == 'arrow':
        run_all_btn['bg'] = '#1e5631'

def on_enter_cancel(e):
    if cancel_btn['cursor'] == 'hand2':
        cancel_btn['bg'] = '#B71C1C'

def on_leave_cancel(e):
    if cancel_btn['cursor'] == 'hand2':
        cancel_btn['bg'] = '#D91E18'

run_all_btn.bind("<Enter>", on_enter_run)
run_all_btn.bind("<Leave>", on_leave_run)
cancel_btn.bind("<Enter>", on_enter_cancel)
cancel_btn.bind("<Leave>", on_leave_cancel)


# ============= FUNCTIONS =============

def run_all():
    global cancel_requested
    cancel_requested = False
    
    filename = file_var.get()
    signatures = sig_var.get()
    sequence = sequence_var.get()
    matrix = matrix_var.get()
    database = database_var.get()
    e_value = eval_var.get()
    num_hits = hits_var.get()
    organism = organism_var.get()
    
    # Extract taxonomy ID if available
    taxonomy_id = None
    if organism and "(taxid:" in organism:
        # Extract the taxonomy ID number
        taxonomy_id = organism.split("(taxid:")[1].split(")")[0].strip()
        # Extract just organism name for BLAST query
        organism = organism.split(" (taxid:")[0].strip()
    
    # Validation
    if not filename:
        messagebox.showerror("Error", "Please enter a file name")
        return
    if not signatures:
        messagebox.showerror("Error", "Please enter signature sequences")
        return
    if not sequence:
        messagebox.showerror("Error", "Please enter a sequence or accession number")
        return
    if not matrix:
        messagebox.showerror("Error", "Please select a matrix")
        return
    if not database:
        messagebox.showerror("Error", "Please select a database")
        return
    try:
        float(e_value)
    except:
        messagebox.showerror("Error", "E-value must be a number")
        return
    try:
        int(num_hits)
    except:
        messagebox.showerror("Error", "Number of hits must be an integer")
        return
    
    run_all_btn.config(bg='#95a5a6', fg='#cccccc', cursor='')
    cancel_btn.config(bg='#D91E18', fg='white', cursor='hand2')
    
    progress_window = tk.Toplevel(window)
    progress_window.title("Processing NEV-BLAST Analysis...")
    progress_window.geometry('700x350')
    progress_window.resizable(False, False)
    progress_window.configure(bg=UWP_WHITE)
    progress_window.transient(window)
    progress_window.grab_set()
    
    header_frame = tk.Frame(progress_window, bg=UWP_GREEN)
    header_frame.pack(fill='x')
    tk.Label(header_frame, text="Processing NEV-BLAST Analysis", 
             font=('Arial', 20, 'bold'), bg=UWP_GREEN, fg=UWP_WHITE).pack(pady=18)
    
    status_label = tk.Label(progress_window, text="Initializing...", 
                           font=('Arial', 16, 'bold'), bg=UWP_WHITE, fg=UWP_GREEN, pady=18)
    status_label.pack()
    
    style = ttk.Style()
    style.configure("UWP.Horizontal.TProgressbar", 
                   background=UWP_GREEN, troughcolor='#ecf0f1', thickness=28)
    
    progress = ttk.Progressbar(progress_window, length=600, mode='determinate',
                              style="UWP.Horizontal.TProgressbar")
    progress.pack(pady=18)
    
    percent_label = tk.Label(progress_window, text="0%", font=('Arial', 16, 'bold'), 
                            bg=UWP_WHITE, fg=UWP_GREEN)
    percent_label.pack()
    
    time_label = tk.Label(progress_window, text="Estimated time: Calculating...", 
                         font=('Arial', 12, 'italic'), bg=UWP_WHITE, fg=UWP_BLACK)
    time_label.pack(pady=12)
    
    detail_label = tk.Label(progress_window, text="", 
                           font=('Arial', 12), bg=UWP_WHITE, fg=UWP_BLACK, wraplength=650)
    detail_label.pack(pady=15)
    
    def update_progress(value, status, detail="", time_est=""):
        progress['value'] = value
        percent_label.config(text=f"{int(value)}%")
        status_label.config(text=status)
        if detail:
            detail_label.config(text=detail)
        if time_est:
            time_label.config(text=f"⏱ Estimated time remaining: {time_est}")
        progress_window.update()
    
    def run_workflow():
        nonlocal organism, taxonomy_id
        start_time = time.time()
        
        try:
            if cancel_requested:
                raise Exception("Cancelled by user")
            
            # Create output subfolder
            import os
            output_dir = os.path.join(os.getcwd(), filename)
            os.makedirs(output_dir, exist_ok=True)
            print(f"[OUTPUT] Created directory: {output_dir}")
                
            update_progress(5, "Step 1/3: Running BLAST Search",
                          f"Searching {database} database...", "10-15 minutes")
            
            time.sleep(1)
            update_progress(10, "Step 1/3: Running BLAST Search", 
                          "Connecting to NCBI servers...", "8-12 minutes")
            
            import ssl
            ssl._create_default_https_context = ssl._create_unverified_context
            from Bio.Blast import NCBIWWW
            
            update_progress(15, "Step 1/3: Running BLAST Search", 
                          f"Querying for {num_hits} matches...", "5-10 minutes")
            
            # Run BLAST query
            if organism and organism.strip():
                result_handle = NCBIWWW.qblast(
                    program="blastp",
                    database=database,
                    sequence=sequence,
                    expect=float(e_value),
                    hitlist_size=int(num_hits),
                    entrez_query=f'"{organism}" [Organism]'
                )
            else:
                result_handle = NCBIWWW.qblast(
                    program="blastp",
                    database=database,
                    sequence=sequence,
                    expect=float(e_value),
                    hitlist_size=int(num_hits)
                )
            
            if cancel_requested:
                raise Exception("Cancelled by user")
                
            update_progress(50, "Step 1/3: Running BLAST Search", 
                          "Receiving results from NCBI...", "3-5 minutes")
            
            # Save XML to subfolder
            xml_path = os.path.join(output_dir, filename + ".xml")
            with open(xml_path, "w") as out_handle:
                out_handle.write(result_handle.read())
            result_handle.close()
            
            update_progress(60, "Step 1/3: BLAST Complete ✓", 
                          f"Saved to {filename}/{filename}.xml", "2-4 minutes")
            time.sleep(0.5)
            
            if cancel_requested:
                raise Exception("Cancelled by user")
                
            update_progress(65, "Step 2/3: Processing Signatures", 
                          "Parsing BLAST results...", "2-3 minutes")
            
            clean_sigs = signatures.replace(''', "'").replace(''', "'")
            clean_sigs = clean_sigs.replace('"', '"').replace('"', '"')
            
            # Change to output directory for backend processing
            original_dir = os.getcwd()
            os.chdir(output_dir)
            
            nev.blastparser(clean_sigs, filename)
            
            if cancel_requested:
                os.chdir(original_dir)
                raise Exception("Cancelled by user")
                
            update_progress(75, "Step 2/3: Processing Signatures", 
                          "Scoring with BLOSUM62...", "1-2 minutes")
            
            nev.hash(filename)
            
            update_progress(90, "Step 2/3: Signatures Complete ✓", 
                          f"Saved to {filename}/{filename}.csv", "10-20 seconds")
            time.sleep(0.5)
            
            if cancel_requested:
                os.chdir(original_dir)
                raise Exception("Cancelled by user")
                
            update_progress(95, "Step 3/3: Generating Visualization", 
                          "Creating 3D plot...", "5-10 seconds")
            
            # Pass taxonomy_id to Plot.threedplot
            Plot.threedplot(filename, taxonomy_id=taxonomy_id)
            
            # Return to original directory
            os.chdir(original_dir)
            
            update_progress(100, "Analysis Complete! ✓", 
                          "Plot opened in browser.", "")
            
            total_time = time.time() - start_time
            minutes = int(total_time // 60)
            seconds = int(total_time % 60)
            
            messagebox.showinfo("Success!", 
                              f"NEV-BLAST completed in {minutes}m {seconds}s\n\n"
                              f"All files saved in folder: {filename}/\n\n"
                              f"Files created:\n"
                              f"• {filename}.xml\n• {filename}.txt\n"
                              f"• {filename}.csv\n• {filename}_visualization.html")
            
            progress_window.destroy()
            
        except Exception as e:
            progress_window.destroy()
            if "Cancelled by user" in str(e):
                messagebox.showinfo("Cancelled", "Analysis cancelled")
            else:
                messagebox.showerror("Error", f"Error:\n{str(e)}")
        
        finally:
            run_all_btn.config(bg='#1e5631', fg='white', cursor='arrow')
            cancel_btn.config(bg='#EE7674', fg='white', cursor='')
    
    thread = threading.Thread(target=run_workflow, daemon=True)
    thread.start()

def cancel_workflow():
    global cancel_requested
    if messagebox.askyesno("Cancel?", "Cancel the analysis?"):
        cancel_requested = True


window.mainloop()