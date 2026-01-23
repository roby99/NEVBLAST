# Create a 2-D or 3-D plot with substitution pattern coloring and clickable points
import re
import plotly.graph_objects as go
import csv
import pandas as pd
import string
import webbrowser
import os
import json

def threedplot(name, taxonomy_id=None):
    """
    Generate interactive 3D and 2D plots - optimized and user-friendly
    taxonomy_id: Optional NCBI taxonomy ID (just the number) for organism filtering
    """
    
    def removecharacters(signature):
        """Extract clean amino acid codes from signature string"""
        remove = string.punctuation
        remove = remove.replace("-", "")
        pattern = r"[{}]".format(remove)
        signaturetext = re.sub(pattern, "", signature)
        signaturetext = re.sub(r'[0-9\.]+', '', signaturetext)
        signaturetext = signaturetext.replace(" ", "")
        signaturetext = ', '.join(a+b for a,b in zip(signaturetext[::2], signaturetext[1::2]))
        return signaturetext

    # Data containers
    data_rows = []
    has_3d = False
    
    # Read CSV and extract data
    with open(name + ".csv", "r") as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        r = 0
        for row in plots:
            if r < 2:
                r += 1
                continue
            
            try:
                if len(row) < 6:
                    continue
                
                scores = []
                for i in range(len(row)-1, 2, -1):
                    try:
                        val = float(row[i])
                        scores.insert(0, val)
                        if len(scores) == 2:
                            break
                    except ValueError:
                        continue
                
                sig1_text = removecharacters(row[3]) if len(row) > 3 else ""
                sig2_text = removecharacters(row[4]) if len(row) > 4 else ""
                substitution_pattern = f"{sig1_text} | {sig2_text}" if sig2_text else sig1_text
                
                if len(scores) >= 2:
                    data_rows.append({
                        'sig1_score': scores[0],
                        'sig2_score': scores[1],
                        'evalue': float(row[2]),
                        'accession': row[0],
                        'sequence': row[1] if len(row) > 1 else "",
                        'sig1_text': sig1_text,
                        'sig2_text': sig2_text,
                        'substitution_pattern': substitution_pattern,
                        'plot_type': '3D'
                    })
                elif len(scores) == 1:
                    data_rows.append({
                        'sig1_score': scores[0],
                        'sig2_score': None,
                        'evalue': float(row[2]),
                        'accession': row[0],
                        'sequence': row[1] if len(row) > 1 else "",
                        'sig1_text': sig1_text,
                        'sig2_text': '',
                        'substitution_pattern': substitution_pattern,
                        'plot_type': '2D'
                    })
                    
            except (ValueError, IndexError) as e:
                print(f"Skipping row due to error: {e}")
                continue
    
    if not data_rows:
        print("No valid data found in CSV")
        return
    
    df = pd.DataFrame(data_rows)
    has_3d = any(df['plot_type'] == '3D')
    
    # Create color mapping
    unique_patterns = df['substitution_pattern'].unique()
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
              '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
              '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5']
    
    pattern_colors = {pattern: colors[i % len(colors)] for i, pattern in enumerate(unique_patterns)}
    df['color'] = df['substitution_pattern'].map(pattern_colors)
    
    representative_sequence = df.iloc[0]['sequence'] if len(df) > 0 else ""
    taxonomy_js = str(taxonomy_id) if taxonomy_id else ""
    
    html_parts = []
    
    html_parts.append("""
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>NEV-BLAST Interactive Visualization</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        * {
            box-sizing: border-box;
        }
        html, body {
            margin: 0;
            padding: 0;
            height: 100%;
        }
        
        /* Page */
        body {
            font-family: Arial, sans-serif;
            background-color: #ffffff;
            overflow-y: auto;
        }
        .header {
            background-color: #1e5631;
            color: white;
            padding: 12px;
            text-align: center;
        }
        .header h1 {
            margin: 0 0 5px 0;
            font-size: 22px;
        }
        .header p {
            margin: 0;
            font-size: 18px;
            font-weight: bold;
        }
        .header .selected-text {
            color: #ffffff;
            font-size: 18px;
            font-weight: bold;
            font-style: italic;
            text-decoration: underline;
        }
        .controls {
            background-color: white;
            padding: 12px;
            text-align: center;
            border-bottom: 2px solid #1e5631;
            display: flex;
            justify-content: center;
            align-items: center;
            gap: 10px;
            flex-wrap: wrap;
        }
        .control-group {
            display: flex;
            align-items: center;
            gap: 6px;
        }
        .controls label {
            font-size: 14px;
            font-weight: bold;
        }
        .controls select {
            font-size: 14px;
            padding: 7px 12px;
            border: 1px solid #ccc;
            border-radius: 4px;
            transition: all 0.2s ease;
        }
        .controls button {
            font-size: 14px;
            padding: 7px 14px;
            border: 1px solid #ccc;
            border-radius: 4px;
            transition: all 0.2s ease;
            background-color: #1e5631;
            color: white;
            cursor: pointer;
            border: none;
            font-weight: bold;
        }
        .controls button:hover {
            background-color: #2d7a4a;
            transform: translateY(-1px);
        }
        .controls button:active {
            transform: translateY(0);
        }
        .reset-btn {
            background-color: #d62728 !important;
        }
        .reset-btn:hover {
            background-color: #a01f20 !important;
        }
        .plot-controls {
            display: flex;
            gap: 5px;
        }
        .plot-btn {
            padding: 7px 14px;
            font-size: 13px;
        }
        .database-group {
            display: flex;
            gap: 5px;
            flex-wrap: wrap;
        }
        .db-btn {
            padding: 7px 14px;
            font-size: 13px;
            background-color: #2563eb !important;
        }
        .db-btn:hover {
            background-color: #1d4ed8 !important;
        }
        #plotContainer {
            position: relative;
            float: left;
            width: calc(100% - 320px);
            min-height: calc(100vh - 110px);
            background-color: #ffffff;
            overflow: hidden;
        }
        
        #plot {
            width: 100%;
            height: calc(100vh - 180px);
            margin-top: -30px;
            position: relative;
        }
        
        .loading-overlay {
            position: absolute;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            background-color: rgba(255, 255, 255, 0.8);
            display: none;
            justify-content: center;
            align-items: center;
            z-index: 1000;
        }
        .loading-overlay.active {
            display: flex;
        }
        .spinner {
            border: 4px solid #f3f3f3;
            border-top: 4px solid #1e5631;
            border-radius: 50%;
            width: 40px;
            height: 40px;
            animation: spin 1s linear infinite;
        }
        @keyframes spin {
            0% { transform: rotate(0deg); }
            100% { transform: rotate(360deg); }
        }
        
        /* Fixed panel - top will be set by JavaScript */
        .filter-panel {
            position: fixed;
            top: 0;  /* will be dynamically set by positionFilterPanel() */
            right: 0;
            bottom: 0;
            width: 320px;
            background-color: #20552a;
            border-left: 2px solid #1e5631;
            border-top: 2px solid #1e5631; /* connects to selection bar line */
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
            z-index: 100;
            display: flex;
            flex-direction: column;
            margin: 0;  /* ensure no margin */
        }
        
        .filter-panel h3 {
            margin: 0;
            padding: 2px 10px 6px 10px;  /* minimal top padding for flush appearance */
            font-size: 22px;
            font-weight: 600;
            color: #ffffff;
            border-bottom: 2px solid #ffffff;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }
        
        #font-controls {
            display: flex;
            gap: 4px;
            align-items: center;
        }
        
        #font-controls button {
            background-color: #ffffff;
            color: #20552a;
            border: none;
            border-radius: 3px;
            padding: 2px 6px;
            cursor: pointer;
            font-size: 11px;
            font-weight: 600;
            transition: background-color 0.2s;
        }
        
        #font-controls button:hover {
            background-color: #e5e5e5;
        }
        
        #filterOptions {
            flex: 1;
            overflow-y: auto;
            padding: 6px;
            background-color: #ffffff;
            color: #000000;
            font-size: 14px;
        }
        
        .filter-actions {
            padding: 6px;
            border-top: 1px solid #ccc;
            background-color: #ffffff;
            display: flex;
            gap: 5px;
        }
        
        body::-webkit-scrollbar,
        #filterOptions::-webkit-scrollbar {
            width: 14px;
            height: 14px;
        }
        
        body::-webkit-scrollbar-track,
        #filterOptions::-webkit-scrollbar-track {
            background: #dcdcdc;
            border-radius: 7px;
        }
        
        body::-webkit-scrollbar-thumb,
        #filterOptions::-webkit-scrollbar-thumb {
            background-color: #555;
            border-radius: 7px;
            border: 2px solid #dcdcdc;
        }
        
        body::-webkit-scrollbar-thumb:hover,
        #filterOptions::-webkit-scrollbar-thumb:hover {
            background-color: #333;
        }
        
        body::-webkit-scrollbar-thumb:active,
        #filterOptions::-webkit-scrollbar-thumb:active {
            background-color: #1e5631;
        }
        
        #filterOptions .filter-option label {
            font-size: inherit;
            line-height: 1.4;
        }
        
        .filter-option {
            display: flex;
            align-items: center;
            gap: 8px;
            margin: 6px 0;
            cursor: pointer;
        }
        .filter-option input[type="checkbox"] {
            width: 16px;
            height: 16px;
            cursor: pointer;
            flex-shrink: 0;
        }
        .filter-option label {
            font-size: 13px;
            cursor: pointer;
            flex: 1;
            word-wrap: break-word;
            line-height: 1.3;
        }
        .color-marker {
            width: 14px;
            height: 14px;
            border-radius: 50%;
            display: inline-block;
            margin-right: 5px;
            flex-shrink: 0;
        }
        .filter-actions {
            margin: 10px 6px;
            padding-top: 10px;
            border-top: 1px solid #ccc;
            display: flex;
            gap: 5px;
            flex-shrink: 0;
            background-color: #ffffff;
        }
        .filter-actions button {
            flex: 1;
            padding: 5px;
            font-size: 11px;
        }
        .info-icon {
            display: inline-block;
            width: 16px;
            height: 16px;
            background-color: #1e5631;
            color: white;
            border-radius: 50%;
            text-align: center;
            line-height: 16px;
            font-size: 12px;
            cursor: help;
            margin-left: 5px;
        }
        .tooltip {
            position: relative;
            display: inline-block;
        }
        .tooltiptext {
            visibility: hidden;
            width: 300px;
            background-color: #555;
            color: #fff;
            text-align: left;
            border-radius: 6px;
            padding: 10px;
            position: absolute;
            z-index: 1;
            bottom: 125%;
            left: 50%;
            margin-left: -150px;
            opacity: 0;
            transition: opacity 0.3s;
            font-size: 12px;
            line-height: 1.4;
        }
        .tooltip:hover .tooltiptext {
            visibility: visible;
            opacity: 1;
        }
        
        /* Global Font Size Increase */
        body {
            font-size: 19px;
        }
        
        h1, h2, h3, h4, h5 {
            font-size: 1.3em;
        }
        
        button, select, label, .tooltiptext, .filter-option label, .header p, .header h1 {
            font-size: 1.15em;
        }
        
        .tooltiptext {
            font-size: 1.05em;
            line-height: 1.5;
        }
        
        .filter-panel h3 {
            font-size: 1.25em;
        }
        
        .filter-option {
            font-size: 1.1em;
        }
        
        .plot-btn, .db-btn {
            font-size: 1.1em;
            padding: 10px 18px;
        }
        
        .controls select, .controls label {
            font-size: 1.1em;
        }
        
        .header, .controls {
            padding: 14px;
        }
        
        .selection-bar {
            display: flex;
            flex-wrap: wrap;
            align-items: center;
            gap: 12px;
            padding: 6px 16px;
            padding-right: 336px; /* 320px panel + a bit of spacing */
            border-top: 2px solid #1e5631;
            border-bottom: 2px solid #1e5631;
            background-color: #ffffff;
            font-size: 14px;
            width: 100%;
            box-sizing: border-box;
        }
        
        .sel-label {
            font-weight: 700;
            color: #1e5631;
        }
        
        .sel-subs {
            max-width: 320px;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
            cursor: default;
        }
        
        @media (max-width: 1000px) {
            #plotContainer { width: 100%; }
            .filter-panel { position: relative; width: 100%; border-left: none; box-shadow: none; }
            .selection-info { margin: 8px 20px; }
        }
    </style>
</head>
<body>
    <div class="header">
        <h1>NEV-BLAST Interactive Visualization</h1>
    </div>
    
    <div class="controls">
        <div class="control-group">
            <label for="viewSelector">View:</label>
            <select id="viewSelector" onchange="changeView()">
""")
    
    if has_3d:
        html_parts.append('                <option value="3d">3D Plot (E-value vs Both Signatures)</option>\n')
        html_parts.append('                <option value="2d_sig1">2D Plot (E-value vs Signature 1)</option>\n')
        html_parts.append('                <option value="2d_sig2">2D Plot (E-value vs Signature 2)</option>\n')
    else:
        html_parts.append('                <option value="2d_sig1">2D Plot (E-value vs Signature 1)</option>\n')
    
    html_parts.append("""
            </select>
            
            <button class="reset-btn" onclick="resetView()" title="Reset everything">Reset View</button>
            
            <label>Tools:</label>
            <button class="plot-btn" onclick="zoomIn()">Zoom In</button>
            <button class="plot-btn" onclick="zoomOut()">Zoom Out</button>
            <button class="plot-btn" onclick="downloadPlot()">Download</button>
            
            <label>Protein:</label>
            <button class="db-btn" onclick="openNCBI()">NCBI</button>
            <button class="db-btn" onclick="openInterPro()">InterPro</button>
        </div>
    </div>
    
    <div class="selection-bar">
        <span class="sel-label">Selected:</span>
        <span id="selAccession">None</span>
        
        <span class="sel-label">E-value:</span>
        <span id="selEvalue">—</span>
        
        <span class="sel-label">Sig1:</span>
        <span id="selSig1">—</span>
        
        <span class="sel-label">Sig2:</span>
        <span id="selSig2">—</span>
        
        <span class="sel-label">Subs:</span>
        <span id="selSubstitutions" class="sel-subs" title="—">—</span>
    </div>
    
    <div class="filter-panel" id="filterPanel">
        <h3>
            Substitutions
            <span id="font-controls">
                <button onclick="adjustFontSize(-1)">Font−</button>
                <button onclick="adjustFontSize(1)">Font+</button>
            </span>
        </h3>
        
        <div id="filterOptions"></div>
        
        <div class="filter-actions">
            <button onclick="selectAllFilters()">Select All</button>
            <button onclick="deselectAllFilters()">Clear All</button>
        </div>
    </div>
    
    <div id="plotContainer">
        <div id="plot"></div>
        <div class="loading-overlay" id="loadingOverlay">
            <div class="spinner"></div>
        </div>
    </div>
    
    <script>
        let selectedAccession = null;
        let selectedSequence = """ + json.dumps(representative_sequence) + """;
        let originalTraceData = [];
        let baseColors = [];
        let isUpdating = false;
        let selectedFilters = new Set();
        let currentView = null;
        let fixedLayouts = {
            '3d': null,
            '2d_sig1': null,
            '2d_sig2': null
        };
        
        const SELECTED_MARKER_SIZE_2D = 24;
        const SELECTED_MARKER_SIZE_3D = 14;
        
        const taxonomyId = '""" + taxonomy_js + """';
        
        function showLoading() {
            document.getElementById('loadingOverlay').classList.add('active');
        }
        
        function hideLoading() {
            document.getElementById('loadingOverlay').classList.remove('active');
        }
        
        function updateSelectionInfo(details) {
            if (!details) {
                document.getElementById('selAccession').textContent = 'None';
                document.getElementById('selSubstitutions').textContent = '—';
                document.getElementById('selSubstitutions').title = '—';
                document.getElementById('selEvalue').textContent = '—';
                document.getElementById('selSig1').textContent = '—';
                document.getElementById('selSig2').textContent = '—';
                return;
            }
            
            document.getElementById('selAccession').textContent = details.accession || '—';
            const subsEl = document.getElementById('selSubstitutions');
            subsEl.textContent = details.subs || '—';
            subsEl.title = details.subs || '—';
            document.getElementById('selEvalue').textContent = details.evalue || '—';
            document.getElementById('selSig1').textContent = 
                (details.sig1 !== undefined && details.sig1 !== null) ? details.sig1 : '—';
            document.getElementById('selSig2').textContent = 
                (details.sig2 !== undefined && details.sig2 !== null) ? details.sig2 : '—';
        }
        
        function getGraphDiv() {
            return document.getElementById('plot');
        }
        
        function positionFilterPanel() {
            const selBar = document.querySelector('.selection-bar');
            const panel = document.getElementById('filterPanel');
            if (!selBar || !panel) return;
            
            // getBoundingClientRect().top = distance from top of viewport
            const rect = selBar.getBoundingClientRect();
            panel.style.top = rect.top + 'px';
        }
        
        function zoom(factor) {
            const gd = getGraphDiv();
            if (!gd || !gd._fullLayout) return;
            
            if (currentView === '3d' && gd._fullLayout.scene && gd._fullLayout.scene.camera) {
                const c = gd._fullLayout.scene.camera;
                const newEye = {
                    x: c.eye.x * factor,
                    y: c.eye.y * factor,
                    z: c.eye.z * factor
                };
                Plotly.relayout(gd, { 'scene.camera.eye': newEye });
                return;
            }
            
            const xa = gd._fullLayout.xaxis;
            const ya = gd._fullLayout.yaxis;
            if (!xa || !ya || !xa.range || !ya.range) return;
            
            const xMid = (xa.range[0] + xa.range[1]) / 2;
            const yMid = (ya.range[0] + ya.range[1]) / 2;
            
            const xHalf = (xa.range[1] - xa.range[0]) * factor / 2;
            const yHalf = (ya.range[1] - ya.range[0]) * factor / 2;
            
            Plotly.relayout(gd, {
                'xaxis.autorange': false,
                'yaxis.autorange': false,
                'xaxis.range': [xMid - xHalf, xMid + xHalf],
                'yaxis.range': [yMid - yHalf, yMid + yHalf]
            });
        }
        
        function zoomIn() {
            zoom(0.8);
        }
        
        function zoomOut() {
            zoom(1.25);
        }
        
        function downloadPlot() {
            const gd = getGraphDiv();
            if (!gd) return;
            Plotly.downloadImage(gd, {
                format: 'png',
                width: 1920,
                height: 1080,
                filename: 'NEV-BLAST_plot'
            });
        }
        
        function resetView() {
            changeView();
        }
        
        function openNCBI() {
            if (!selectedAccession) {
                alert('Click a point on the plot first to select a protein');
                return;
            }
            window.open(`https://www.ncbi.nlm.nih.gov/protein/${selectedAccession}`, '_blank');
        }
        
        function openInterPro() {
            if (!selectedSequence) {
                alert('Click a point on the plot first to select a protein');
                return;
            }
            window.open(`https://www.ebi.ac.uk/interpro/search/sequence/${encodeURIComponent(selectedSequence)}`, '_blank');
        }
        
        function createFilterPanel() {
            const filterOptions = document.getElementById('filterOptions');
            filterOptions.innerHTML = '';
            selectedFilters = new Set();
            
            for (let i = 0; i < originalTraceData.length; i++) {
                const trace = originalTraceData[i];
                selectedFilters.add(i);
                
                const option = document.createElement('div');
                option.className = 'filter-option';
                
                const checkbox = document.createElement('input');
                checkbox.type = 'checkbox';
                checkbox.id = `filter_${i}`;
                checkbox.checked = true;
                checkbox.onchange = () => toggleFilter(i);
                
                const label = document.createElement('label');
                label.htmlFor = `filter_${i}`;
                label.innerHTML = `<span class="color-marker" style="background-color: ${baseColors[i]};"></span>${trace.name}`;
                
                option.appendChild(checkbox);
                option.appendChild(label);
                filterOptions.appendChild(option);
            }
        }
        
        function adjustFontSize(change) {
            const list = document.getElementById('filterOptions');
            const style = window.getComputedStyle(list, null).getPropertyValue('font-size');
            const currentSize = parseFloat(style);
            list.style.fontSize = (currentSize + change) + 'px';
        }
        
        function toggleFilter(traceIndex) {
            const checkbox = document.getElementById(`filter_${traceIndex}`);
            if (checkbox.checked) {
                selectedFilters.add(traceIndex);
            } else {
                selectedFilters.delete(traceIndex);
            }
            applyFilters();
        }
        
        function selectAllFilters() {
            selectedFilters.clear();
            for (let i = 0; i < originalTraceData.length; i++) {
                selectedFilters.add(i);
                document.getElementById(`filter_${i}`).checked = true;
            }
            applyFilters();
        }
        
        function deselectAllFilters() {
            selectedFilters.clear();
            for (let i = 0; i < originalTraceData.length; i++) {
                document.getElementById(`filter_${i}`).checked = false;
            }
            applyFilters();
        }
        
        function applyFilters() {
            if (isUpdating) return;
            isUpdating = true;
            
            const gd = getGraphDiv();
            const newData = [];
            
            for (let i = 0; i < originalTraceData.length; i++) {
                if (!selectedFilters.has(i)) continue;
                newData.push(originalTraceData[i]);
            }
            
            if (newData.length === 0 && originalTraceData.length > 0) {
                newData.push(originalTraceData[0]);
            }
            
            let layout;
            if (currentView === '3d') {
                layout = layout3d;
            } else if (currentView === '2d_sig2') {
                layout = layout2d_sig2;
            } else {
                layout = layout2d_sig1;
            }
            
            Plotly.newPlot(gd, newData, layout, plotConfig)
                .then(() => {
                    gd.on('plotly_click', onPointClick);
                    gd.on('plotly_legendclick', onLegendClick);
                    gd.on('plotly_legenddoubleclick', onLegendDoubleClick);
                    isUpdating = false;
                })
                .catch(() => {
                    isUpdating = false;
                });
        }
        
        function captureFixedLayout() {
            const gd = getGraphDiv();
            const full = gd._fullLayout;
            if (!full || !currentView) return;
            
            if (currentView === '3d' && full.scene) {
                const copy = JSON.parse(JSON.stringify(layout3d));
                copy.scene.domain = full.scene.domain;
                copy.scene.xaxis.range = full.scene.xaxis.range.slice();
                copy.scene.yaxis.range = full.scene.yaxis.range.slice();
                copy.scene.zaxis.range = full.scene.zaxis.range.slice();
                copy.scene.xaxis.autorange = false;
                copy.scene.yaxis.autorange = false;
                copy.scene.zaxis.autorange = false;
                fixedLayouts['3d'] = copy;
            } else if (currentView === '2d_sig1' && full.xaxis && full.yaxis) {
                const copy = JSON.parse(JSON.stringify(layout2d_sig1));
                copy.xaxis.range = full.xaxis.range.slice();
                copy.yaxis.range = full.yaxis.range.slice();
                copy.xaxis.autorange = false;
                copy.yaxis.autorange = false;
                fixedLayouts['2d_sig1'] = copy;
            } else if (currentView === '2d_sig2' && full.xaxis && full.yaxis) {
                const copy = JSON.parse(JSON.stringify(layout2d_sig2));
                copy.xaxis.range = full.xaxis.range.slice();
                copy.yaxis.range = full.yaxis.range.slice();
                copy.xaxis.autorange = false;
                copy.yaxis.autorange = false;
                fixedLayouts['2d_sig2'] = copy;
            }
        }
        
        function onPointClick(data) {
            if (isUpdating || !data.points || data.points.length === 0) return;
            
            const point = data.points[0];
            if (!point.customdata || point.customdata.length < 2) return;
            
            const accession = point.customdata[0];
            const sequence = point.customdata[1];
            
            const gd = getGraphDiv();
            if (gd.data &&
                gd.data.length === 1 &&
                gd.data[0].customdata &&
                gd.data[0].customdata[0] &&
                gd.data[0].customdata[0][0] === accession) {
                changeView();
                return;
            }
            
            isUpdating = true;
            
            selectedAccession = accession;
            selectedSequence = sequence;
            
            const baseTrace = point.data;
            
            let subs = '';
            let evalue = '';
            let sig1 = '';
            let sig2 = '';
            
            if (baseTrace && baseTrace.name) {
                subs = baseTrace.name;
            }
            if (point.y !== undefined) {
                evalue = Array.isArray(point.y) ? point.y[0] : point.y;
            }
            if (point.x !== undefined) {
                sig1 = Array.isArray(point.x) ? point.x[0] : point.x;
            }
            if (point.z !== undefined) {
                sig2 = Array.isArray(point.z) ? point.z[0] : point.z;
            }
            
            updateSelectionInfo({
                accession,
                subs,
                evalue: evalue !== '' ? evalue.toExponential ? evalue.toExponential(2) : evalue : '—',
                sig1,
                sig2
            });
            
            const is3D = baseTrace.type === 'scatter3d';
            const singleTrace = {
                type: baseTrace.type,
                mode: 'markers',
                name: baseTrace.name,
                x: [point.x],
                y: [point.y],
                marker: {
                    size: is3D ? SELECTED_MARKER_SIZE_3D : SELECTED_MARKER_SIZE_2D,
                    color: (baseTrace.marker && baseTrace.marker.color)
                        ? baseTrace.marker.color
                        : '#d62728',
                    opacity: 1,
                    line: { width: 4, color: '#222' }
                },
                text: [
                    baseTrace.text
                        ? baseTrace.text[point.pointIndex] || ''
                        : ''
                ],
                hoverinfo: baseTrace.hoverinfo || 'text',
                customdata: [[accession, sequence]]
            };
            
            if (baseTrace.type === 'scatter3d' && point.z !== undefined) {
                singleTrace.z = [point.z];
            }
            
            let layout;
            if (fixedLayouts[currentView]) {
                layout = JSON.parse(JSON.stringify(fixedLayouts[currentView]));
            } else {
                if (currentView === '3d') {
                    layout = JSON.parse(JSON.stringify(layout3d));
                } else if (currentView === '2d_sig2') {
                    layout = JSON.parse(JSON.stringify(layout2d_sig2));
                } else {
                    layout = JSON.parse(JSON.stringify(layout2d_sig1));
                }
            }
            
            if (layout.yaxis && layout.yaxis.type === 'log' && singleTrace.y[0] <= 0) {
                singleTrace.y[0] = 1e-300;
            }
            if (layout.scene && layout.scene.zaxis && layout.scene.zaxis.type === 'log' 
                && singleTrace.z && singleTrace.z[0] <= 0) {
                singleTrace.z[0] = 1e-300;
            }
            
            gd.style.opacity = 0;
            setTimeout(() => {
                Plotly.newPlot(gd, [singleTrace], layout, plotConfig).then(() => {
                    gd.style.opacity = 1;
                    gd.on('plotly_click', onPointClick);
                    gd.on('plotly_legendclick', onLegendClick);
                    gd.on('plotly_legenddoubleclick', onLegendDoubleClick);
                    isUpdating = false;
                });
            }, 150);
        }
        
        function onLegendClick(data) {
            return false;
        }
        
        function onLegendDoubleClick() {
            return false;
        }
        
        function storeOriginalData(data) {
            originalTraceData = JSON.parse(JSON.stringify(data));
            baseColors = [];
            for (let i = 0; i < data.length; i++) {
                baseColors.push(data[i].marker.color);
            }
            
            createFilterPanel();
        }
""")
    
    shared_margins = {'l': 10, 'r': 10, 't': 10, 'b': 10}
    # Plot height removed - now handled dynamically by CSS calc(100vh - 260px)
    
    if has_3d:
        df_3d = df[df['plot_type'] == '3D']
        traces_3d = []
        
        for pattern in unique_patterns:
            df_pattern = df_3d[df_3d['substitution_pattern'] == pattern]
            if len(df_pattern) > 0:
                hover_text = [f"<b>{row['accession']}</b><br>Sig1: {row['sig1_text']}<br>Sig2: {row['sig2_text']}<br>E-val: {row['evalue']:.2e}" 
                             for _, row in df_pattern.iterrows()]
                
                trace = {
                    'type': 'scatter3d',
                    'mode': 'markers',
                    'name': pattern,
                    'x': df_pattern['sig1_score'].tolist(),
                    'y': df_pattern['sig2_score'].tolist(),
                    'z': df_pattern['evalue'].tolist(),
                    'marker': {
                        'size': 5,
                        'color': pattern_colors[pattern],
                        'opacity': 0.85,
                        'line': {'width': 0}
                    },
                    'text': hover_text,
                    'hoverinfo': 'text',
                    'customdata': [[row['accession'], row['sequence']] for _, row in df_pattern.iterrows()]
                }
                traces_3d.append(trace)
        
        layout_3d = {
            'autosize': True,
            'scene': {
                'domain': {
                    'y': [-0.1, .99]
                },
                'xaxis': {
                    'title': 'Signature 1',
                    'titlefont': {'size': 16}
                },
                'yaxis': {
                    'title': 'Signature 2',
                    'titlefont': {'size': 16}
                },
                'zaxis': {
                    'title': 'E-value',
                    'type': 'log',
                    'titlefont': {'size': 16}
                },
                'camera': {
                    'eye': {
                            'x': 1.706917760350791,
                            'y': 2.049526994724092,
                            'z': 0.16277548172966846
                            },
                    'center': {'x': 0, 'y': 0, 'z': 0},
                    'up': {'x': 0, 'y': 0, 'z': 1}
                },
                'aspectmode': 'manual',
                'aspectratio': {'x': 2.5, 'y': 1, 'z': 1}
            },
            'title': {
                'text': '3D Visualization: E-value vs Signature Substitutions',
                'y': 0.94,
                'x': 0.5,
                'xanchor': 'center',
                'yanchor': 'top',
                'font': {'size': 20}
            },
            'showlegend': False,
            'hovermode': 'closest',
            'hoverlabel': {'font': {'size': 15}},
            'margin': shared_margins
        }
        
        html_parts.append(f"        const data3d = {json.dumps(traces_3d)};\n")
        html_parts.append(f"        const layout3d = {json.dumps(layout_3d)};\n")
    
    traces_2d_sig1 = []
    for pattern in unique_patterns:
        df_pattern = df[df['substitution_pattern'] == pattern]
        if len(df_pattern) > 0:
            hover_text = [f"<b>{row['accession']}</b><br>Sig1: {row['sig1_text']}<br>E-val: {row['evalue']:.2e}" 
                         for _, row in df_pattern.iterrows()]
            
            trace = {
                'type': 'scatter',
                'mode': 'markers',
                'name': pattern,
                'x': df_pattern['sig1_score'].tolist(),
                'y': df_pattern['evalue'].tolist(),
                'marker': {
                    'size': 14,
                    'color': pattern_colors[pattern],
                    'opacity': 0.85,
                    'line': {'width': 0}
                },
                'text': hover_text,
                'hoverinfo': 'text',
                'customdata': [[row['accession'], row['sequence']] for _, row in df_pattern.iterrows()]
            }
            traces_2d_sig1.append(trace)
    
    layout_2d_sig1 = {
        'autosize': True,
        'xaxis': {
            'title': 'Signature 1 Score',
            'range': [-0.15, 1.15],
            'titlefont': {'size': 16},
            'tickfont': {'size': 14}
        },
        'yaxis': {
            'title': 'E-value',
            'type': 'log',
            'titlefont': {'size': 16},
            'tickfont': {'size': 14}
        },
        'title': {
            'text': '2D Plot: E-value vs Signature 1',
            'font': {'size': 20}
        },
        'showlegend': False,
        'hovermode': 'closest',
        'hoverlabel': {'font': {'size': 15}},
        'margin': shared_margins
    }
    
    html_parts.append(f"        const data2d_sig1 = {json.dumps(traces_2d_sig1)};\n")
    html_parts.append(f"        const layout2d_sig1 = {json.dumps(layout_2d_sig1)};\n")
    
    if has_3d:
        traces_2d_sig2 = []
        for pattern in unique_patterns:
            df_pattern = df[(df['substitution_pattern'] == pattern) & (df['plot_type'] == '3D')]
            if len(df_pattern) > 0:
                hover_text = [f"<b>{row['accession']}</b><br>Sig2: {row['sig2_text']}<br>E-val: {row['evalue']:.2e}" 
                             for _, row in df_pattern.iterrows()]
                
                trace = {
                    'type': 'scatter',
                    'mode': 'markers',
                    'name': pattern,
                    'x': df_pattern['sig2_score'].tolist(),
                    'y': df_pattern['evalue'].tolist(),
                    'marker': {
                        'size': 14,
                        'color': pattern_colors[pattern],
                        'opacity': 0.85,
                        'line': {'width': 0}
                    },
                    'text': hover_text,
                    'hoverinfo': 'text',
                    'customdata': [[row['accession'], row['sequence']] for _, row in df_pattern.iterrows()]
                }
                traces_2d_sig2.append(trace)
        
        layout_2d_sig2 = {
            'autosize': True,
            'xaxis': {
                'title': 'Signature 2 Score',
                'range': [-0.15, 1.15],
                'titlefont': {'size': 16},
                'tickfont': {'size': 14}
            },
            'yaxis': {
                'title': 'E-value',
                'type': 'log',
                'titlefont': {'size': 16},
                'tickfont': {'size': 14}
            },
            'title': {
                'text': '2D Plot: E-value vs Signature 2',
                'font': {'size': 20}
            },
            'showlegend': False,
            'hovermode': 'closest',
            'hoverlabel': {'font': {'size': 15}},
            'margin': shared_margins
        }
        
        html_parts.append(f"        const data2d_sig2 = {json.dumps(traces_2d_sig2)};\n")
        html_parts.append(f"        const layout2d_sig2 = {json.dumps(layout_2d_sig2)};\n")
    
    html_parts.append("""
        const plotConfig = {
            displayModeBar: false,
            responsive: true
        };
        
        function changeView() {
            showLoading();
            const view = document.getElementById('viewSelector').value;
            const gd = getGraphDiv();
            
            selectedAccession = null;
            updateSelectionInfo(null);
            
            setTimeout(() => {
                if (view === '3d' && typeof data3d !== 'undefined') {
                    currentView = '3d';
                    Plotly.newPlot(gd, data3d, layout3d, plotConfig)
                        .then(() => {
                            storeOriginalData(data3d);
                            captureFixedLayout();
                        });
                } else if (view === '2d_sig1') {
                    currentView = '2d_sig1';
                    Plotly.newPlot(gd, data2d_sig1, layout2d_sig1, plotConfig)
                        .then(() => {
                            storeOriginalData(data2d_sig1);
                            captureFixedLayout();
                        });
                } else if (view === '2d_sig2' && typeof data2d_sig2 !== 'undefined') {
                    currentView = '2d_sig2';
                    Plotly.newPlot(gd, data2d_sig2, layout2d_sig2, plotConfig)
                        .then(() => {
                            storeOriginalData(data2d_sig2);
                            captureFixedLayout();
                        });
                }
                
                gd.on('plotly_click', onPointClick);
                gd.on('plotly_legendclick', onLegendClick);
                gd.on('plotly_legenddoubleclick', onLegendDoubleClick);
                
                positionFilterPanel();
                hideLoading();
            }, 10);
        }
        
        window.addEventListener('load', () => {
            showLoading();
            updateSelectionInfo(null);
            
            setTimeout(() => {
                const gd = getGraphDiv();
                
                if (typeof data3d !== 'undefined') {
                    currentView = '3d';
                    Plotly.newPlot(gd, data3d, layout3d, plotConfig)
                        .then(() => {
                            storeOriginalData(data3d);
                            captureFixedLayout();
                        });
                } else {
                    currentView = '2d_sig1';
                    Plotly.newPlot(gd, data2d_sig1, layout2d_sig1, plotConfig)
                        .then(() => {
                            storeOriginalData(data2d_sig1);
                            captureFixedLayout();
                        });
                }
                
                gd.on('plotly_click', onPointClick);
                gd.on('plotly_legendclick', onLegendClick);
                gd.on('plotly_legenddoubleclick', onLegendDoubleClick);
                
                positionFilterPanel();
                hideLoading();
            }, 10);
        });
        
        window.addEventListener('resize', positionFilterPanel);
        // No scroll listener needed - fixed position stays put while scrolling
    </script>
</body>
</html>
""")
    
    html_content = ''.join(html_parts)
    html_filename = name + "_visualization.html"
    
    with open(html_filename, 'w') as f:
        f.write(html_content)
    
    print(f"Visualization saved to: {html_filename}")
    webbrowser.open('file://' + os.path.abspath(html_filename))