import React, { useState, useEffect, useMemo, useRef } from 'react';
import { Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, Area, ComposedChart, ReferenceArea, ReferenceLine } from 'recharts';
import { Search, Settings, Activity, FileText, Database, Info, AlertCircle, ChevronRight, Sliders, Download, Play, Loader, CheckCircle, Save, Upload } from 'lucide-react';

// --- LOGIC & UTILS ---

// Density Algorithm with Smoothing (Reusable)
const calculateDensity = (positions, length, windowSize, smoothingFactor = 0) => {
  if (!length || length <= 0) return [];
  
  const rawDensity = new Array(length).fill(0);
  const halfWindow = Math.floor(windowSize / 2);

  for (let i = 0; i < length; i++) {
    const start = Math.max(0, i - halfWindow);
    const end = Math.min(length, i + halfWindow);
    let count = 0;
    for (let p of positions) {
      if (p >= start && p <= end) count++;
    }
    rawDensity[i] = count;
  }

  if (smoothingFactor > 0) {
    const smoothedDensity = new Array(length).fill(0);
    const smoothWindow = Math.floor(smoothingFactor);
    
    for (let i = 0; i < length; i++) {
      let sum = 0;
      let count = 0;
      for (let j = i - smoothWindow; j <= i + smoothWindow; j++) {
        if (j >= 0 && j < length) {
          sum += rawDensity[j];
          count++;
        }
      }
      smoothedDensity[i] = sum / count;
    }
    return smoothedDensity;
  }

  return rawDensity;
};

// Pure Analytic Function for a single protein
const analyzeProtein = (protein, params) => {
    // Note: Phosphorylation uses the SAME parameters as Glycosylation (glyco*)
    const { 
        ssWindowSize, ssThreshold, ssSmoothing, 
        glycoWindowSize, glycoThreshold, glycoSmoothing 
    } = params;
    const length = protein.length;

    // 1. Calculate Densities
    const ssDensity = calculateDensity(protein.ssBonds, length, ssWindowSize, ssSmoothing);
    const nDensity = calculateDensity(protein.nLinked, length, glycoWindowSize, glycoSmoothing);
    // Phosphorylation shares glyco parameters
    const pDensity = calculateDensity(protein.phosphorylation, length, glycoWindowSize, glycoSmoothing);
    
    // We define the "High Density Zone" based on SS + N-Linked (the primary analysis target)
    const isHighDensityZone = new Array(length).fill(false);
    for(let i=0; i<length; i++) {
        if(ssDensity[i] >= ssThreshold && nDensity[i] >= glycoThreshold) {
            isHighDensityZone[i] = true;
        }
    }

    // Helper to process a list of sites
    const processSites = (sites, typeCode, typeLabel) => {
        return sites.map(pos => {
            const arrayIndex = pos - 1; 
            
            // Check if embedded in a Disulfide (UniProt Range)
            let embeddedInSS = false;
            if (protein.ssBondRanges) {
                embeddedInSS = protein.ssBondRanges.some(range => pos >= range.start && pos <= range.end);
            }

            // Check if in Co-localization zone
            const inHighDensity = isHighDensityZone[arrayIndex] || false;

            return {
                gene: protein.gene,
                accession: protein.id,
                scientificName: protein.scientificName,
                glycosite: `${typeCode}${pos}`, // e.g. N123, P456
                type: typeLabel,
                inUniProtDisulfide: embeddedInSS,
                inHighDensityZone: inHighDensity
            };
        });
    };

    const nRows = processSites(protein.nLinked, 'N', 'N-Linked');
    const oRows = processSites(protein.oLinked, 'O', 'O-Linked');
    const pRows = processSites(protein.phosphorylation, 'P', 'Phospho');

    return [...nRows, ...oRows, ...pRows];
};

const fetchUniprotData = async (geneName) => {
  const maxRetries = 2;
  let attempt = 0;

  while (attempt <= maxRetries) {
    try {
      const cleanGene = geneName.trim();
      const query = `(gene_exact:"${cleanGene}") AND (organism_id:9606) AND (reviewed:true)`;
      // FIXED: Using 'ft_mod_res' which is the correct field name for Modified Residues in UniProt API
      const fields = "accession,gene_names,length,sequence,ft_disulfid,ft_carbohyd,ft_mod_res,organism_name";
      const searchUrl = `https://rest.uniprot.org/uniprotkb/search?query=${encodeURIComponent(query)}&fields=${fields}&format=json&size=1`;
      
      const searchRes = await fetch(searchUrl);
      
      if (!searchRes.ok) {
        let errorDetails = searchRes.statusText;
        try {
            const errorJson = await searchRes.json();
            if (errorJson.messages) errorDetails = errorJson.messages.join(", ");
        } catch (e) { /* ignore */ }
        throw new Error(`UniProt API Error (${searchRes.status}): ${errorDetails}`);
      }
      
      const searchData = await searchRes.json();

      if (!searchData.results || searchData.results.length === 0) {
        throw new Error(`Gene "${geneName}" not found in Homo sapiens (TaxID: 9606).`);
      }

      return parseUniProtEntry(searchData.results[0]);
    } catch (error) {
      console.error(`Fetch Error (Attempt ${attempt + 1}/${maxRetries + 1}):`, error);
      if (error.message.includes("not found")) return null;
      attempt++;
      if (attempt > maxRetries) return null;
      await new Promise(resolve => setTimeout(resolve, 1000));
    }
  }
  return null;
};

// Centralized parser to convert UniProt JSON to our internal format
const parseUniProtEntry = (entry) => {
    const sequence = entry.sequence.value;
    const length = entry.sequence.length;
    const features = entry.features || [];
    
    const scientificName = entry.organism ? entry.organism.scientificName : "Unknown Species";

    let geneName = entry.primaryAccession;
    if (entry.genes && entry.genes.length > 0 && entry.genes[0].geneName) {
        geneName = entry.genes[0].geneName.value;
    }

    const ssBonds = [];
    const ssBondRanges = [];
    const nLinked = [];
    const oLinked = [];
    const phosphorylation = [];

    features.forEach(f => {
        if (f.type === 'Disulfide bond') {
            ssBonds.push(f.location.start.value);
            ssBondRanges.push({ start: f.location.start.value, end: f.location.end.value });
        } else if (f.type === 'Glycosylation') {
            if (f.description.includes('N-linked')) {
                nLinked.push(f.location.start.value);
            } else if (f.description.includes('O-linked')) {
                oLinked.push(f.location.start.value);
            }
        } else if (f.type === 'Modified residue') {
            const desc = f.description ? f.description.toLowerCase() : "";
            if (desc.includes('phospho')) {
                phosphorylation.push(f.location.start.value);
            }
        }
    });

    return {
        id: entry.primaryAccession,
        gene: geneName,
        scientificName,
        length,
        sequence,
        ssBonds,
        ssBondRanges, 
        nLinked,
        oLinked,
        phosphorylation
    };
};

// Simple FASTA Parser
const parseFasta = (text) => {
    const proteins = [];
    const lines = text.split('\n');
    let currentHeader = null;
    let currentSeq = "";

    const finishEntry = () => {
        if (currentHeader && currentSeq) {
            const nLinked = [];
            const ssBonds = []; 
            const regexN = /N(?=[^P][ST][^P])/g;
            let match;
            while ((match = regexN.exec(currentSeq)) !== null) {
                nLinked.push(match.index + 1); 
            }

            for (let i = 0; i < currentSeq.length; i++) {
                if (currentSeq[i] === 'C') ssBonds.push(i + 1);
            }

            let gene = "Unknown";
            const parts = currentHeader.split('|');
            if (parts.length >= 3) {
                gene = parts[2].split(' ')[0].split('_')[0];
            } else {
                gene = currentHeader.substring(1, 15);
            }

            proteins.push({
                id: currentHeader.substring(1, 10),
                gene: gene,
                scientificName: "Custom FASTA",
                length: currentSeq.length,
                sequence: currentSeq,
                ssBonds: ssBonds, 
                ssBondRanges: [], 
                nLinked: nLinked,
                oLinked: [],
                phosphorylation: []
            });
        }
    };

    for (const line of lines) {
        const l = line.trim();
        if (!l) continue;
        if (l.startsWith('>')) {
            finishEntry();
            currentHeader = l;
            currentSeq = "";
        } else {
            currentSeq += l;
        }
    }
    finishEntry();
    return proteins;
};

// --- COMPONENTS ---

const Card = ({ children, className = "" }) => (
  <div className={`bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden ${className}`}>
    {children}
  </div>
);

const Badge = ({ color, text }) => {
  const colors = {
    blue: "bg-sky-100 text-sky-800 border-sky-200",
    red: "bg-red-100 text-red-800 border-red-200",
    yellow: "bg-yellow-100 text-yellow-800 border-yellow-200",
    purple: "bg-purple-100 text-purple-800 border-purple-200",
    green: "bg-emerald-100 text-emerald-800 border-emerald-200",
    gray: "bg-slate-100 text-slate-600 border-slate-200"
  };
  return (
    <span className={`px-2 py-0.5 rounded text-xs font-medium border ${colors[color] || colors.gray}`}>
      {text}
    </span>
  );
};

const ParameterControl = ({ label, value, onChange, min, max, step, unit, colorClass = "accent-indigo-600" }) => (
  <div className="mb-3">
    <div className="flex justify-between mb-1">
      <label className="text-xs font-medium text-slate-600">{label}</label>
      <span className="text-xs font-bold text-slate-800">{value} {unit}</span>
    </div>
    <input 
      type="range" 
      min={min} 
      max={max} 
      step={step}
      value={value} 
      onChange={(e) => onChange(Number(e.target.value))}
      className={`w-full h-1.5 bg-slate-200 rounded-lg appearance-none cursor-pointer ${colorClass}`}
    />
  </div>
);

export default function App() {
  // State
  const [geneName, setGeneName] = useState("FN1");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [data, setData] = useState(null);
  
  // --- Parameters ---
  // Group 1: SS
  const [ssWindowSize, setSsWindowSize] = useState(50);
  const [ssThreshold, setSsThreshold] = useState(3);
  const [ssSmoothing, setSsSmoothing] = useState(5); 
  
  // Group 2: Glyco (Shared with Phospho)
  const [glycoWindowSize, setGlycoWindowSize] = useState(50);
  const [glycoThreshold, setGlycoThreshold] = useState(3);
  const [glycoSmoothing, setGlycoSmoothing] = useState(5); 
  
  const [view, setView] = useState("calibration"); 

  // Batch State
  const [organismId, setOrganismId] = useState("9606");
  const [batchProgress, setBatchProgress] = useState(0);
  const [batchStatus, setBatchStatus] = useState("idle"); 
  const [batchResults, setBatchResults] = useState([]);
  const [processedCount, setProcessedCount] = useState(0);
  
  useEffect(() => {
    handleSearch();
  }, []);

  const handleSearch = async () => {
    if (!geneName.trim()) return;
    setLoading(true);
    setError(null);
    
    const result = await fetchUniprotData(geneName);
    
    if (result) {
      setData(result);
    } else {
      setError(`Could not retrieve data for gene "${geneName}". Please check the name and try again.`);
    }
    setLoading(false);
  };

  // --- IMPORT / EXPORT PARAMETERS ---
  const exportParameters = () => {
      const params = {
          ssWindowSize, ssThreshold, ssSmoothing,
          glycoWindowSize, glycoThreshold, glycoSmoothing
      };
      const blob = new Blob([JSON.stringify(params, null, 2)], { type: 'text/plain' });
      const url = URL.createObjectURL(blob);
      const link = document.createElement('a');
      link.href = url;
      link.download = 'protein_coloc_params.txt';
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
  };

  const importParameters = (e) => {
      const file = e.target.files[0];
      if (!file) return;
      const reader = new FileReader();
      reader.onload = (evt) => {
          try {
              const params = JSON.parse(evt.target.result);
              if (params.ssWindowSize) setSsWindowSize(params.ssWindowSize);
              if (params.ssThreshold) setSsThreshold(params.ssThreshold);
              if (params.ssSmoothing) setSsSmoothing(params.ssSmoothing);
              if (params.glycoWindowSize) setGlycoWindowSize(params.glycoWindowSize);
              if (params.glycoThreshold) setGlycoThreshold(params.glycoThreshold);
              if (params.glycoSmoothing) setGlycoSmoothing(params.glycoSmoothing);
              alert("Parameters loaded successfully!");
          } catch (err) {
              alert("Error parsing parameter file.");
          }
      };
      reader.readAsText(file);
  };

  // --- EXPORT SINGLE ANALYSIS LOGIC ---
  const handleExport = () => {
    if (!data) return;

    const ssDensity = calculateDensity(data.ssBonds, data.length, ssWindowSize, ssSmoothing);
    const nDensity = calculateDensity(data.nLinked, data.length, glycoWindowSize, glycoSmoothing);
    const oDensity = calculateDensity(data.oLinked, data.length, glycoWindowSize, glycoSmoothing);
    const pDensity = calculateDensity(data.phosphorylation, data.length, glycoWindowSize, glycoSmoothing); // Uses Glyco params

    let csvContent = "Scientific_Name,Position,Amino_Acid,SS_Density,N_Linked_Density,O_Linked_Density,Phospho_Density,Is_HighDensity_Zone\n";

    for (let i = 0; i < data.length; i++) {
        const isOverlap = (ssDensity[i] >= ssThreshold && nDensity[i] >= glycoThreshold);
        const pos = i + 1; 
        const aa = data.sequence[i] || '-';
        
        const row = [
            `"${data.scientificName}"`,
            pos,
            aa,
            ssDensity[i].toFixed(4), 
            nDensity[i].toFixed(4),
            oDensity[i].toFixed(4),
            pDensity[i].toFixed(4),
            isOverlap ? "TRUE" : "FALSE"
        ].join(",");
        
        csvContent += row + "\n";
    }

    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.setAttribute('download', `${geneName}_colocalization_analysis.csv`);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  };


  // --- BATCH PROCESSING LOGIC ---

  const downloadBatchCSV = (rows) => {
      const header = "Scientific_Name,Gene_Name,Type,Site_ID,In_UniProt_Disulfide_Bond,In_HighDensity_Colocalization_Zone\n";
      const content = rows.map(r => 
        `"${r.scientificName}",${r.gene},${r.type},${r.glycosite},${r.inUniProtDisulfide ? "TRUE" : "FALSE"},${r.inHighDensityZone ? "TRUE" : "FALSE"}`
      ).join('\n');
      
      const blob = new Blob([header + content], { type: 'text/csv;charset=utf-8;' });
      const url = URL.createObjectURL(blob);
      const link = document.createElement('a');
      link.href = url;
      link.setAttribute('download', `proteome_analysis_results.csv`);
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
  };

  const processBatchList = (proteins) => {
      // Pass glyco params as they are shared for phospho
      const params = { ssWindowSize, ssThreshold, ssSmoothing, glycoWindowSize, glycoThreshold, glycoSmoothing };
      let results = [];
      proteins.forEach(prot => {
          const analysis = analyzeProtein(prot, params);
          results = [...results, ...analysis];
      });
      return results;
  };

  const runWholeOrganismAnalysis = async () => {
      setBatchStatus("processing");
      setBatchProgress(0);
      setProcessedCount(0);
      setBatchResults([]);
      setError(null);

      let allReportRows = [];
      let nextCursor = null;
      let totalProcessed = 0;
      const size = 250; 
      let hasMore = true;

      try {
          while(hasMore) {
             const query = `(organism_id:${organismId}) AND (reviewed:true)`;
             // FIXED: Using 'ft_mod_res' here as well
             const fields = "accession,gene_names,length,sequence,ft_disulfid,ft_carbohyd,ft_mod_res,organism_name";
             let url = `https://rest.uniprot.org/uniprotkb/search?query=${encodeURIComponent(query)}&fields=${fields}&format=json&size=${size}`;
             if (nextCursor) url += `&cursor=${nextCursor}`;

             const res = await fetch(url);
             if (!res.ok) {
                 let errorMsg = res.statusText;
                 try { const e = await res.json(); if(e.messages) errorMsg = e.messages.join(", "); } catch(x){}
                 throw new Error(`API Error ${res.status}: ${errorMsg}`);
             }
             const json = await res.json();

             const total = Number(res.headers.get('x-total-results')) || 20400;

             const proteins = json.results.map(parseUniProtEntry);
             const chunkResults = processBatchList(proteins);
             allReportRows = [...allReportRows, ...chunkResults];

             totalProcessed += proteins.length;
             setProcessedCount(totalProcessed);
             setBatchProgress(Math.min(100, (totalProcessed / total) * 100));

             if (json.results.length < size || !res.headers.get('link')) {
                 hasMore = false;
             } else {
                 const linkHeader = res.headers.get('link');
                 if (linkHeader && linkHeader.includes('rel="next"')) {
                     const match = linkHeader.match(/[?&]cursor=([^&>"]+)/);
                     if (match) nextCursor = match[1];
                     else hasMore = false;
                 } else {
                     hasMore = false;
                 }
             }
             
             if (totalProcessed > 5000 && organismId === "9606") {
                 // Limit for demo
             }
          }

          setBatchStatus("completed");
          setBatchResults(allReportRows);
          downloadBatchCSV(allReportRows);

      } catch (err) {
          console.error(err);
          setError(`Error during bulk analysis: ${err.message}`);
          setBatchStatus("idle");
      }
  };

  const handleFastaUpload = (e) => {
      const file = e.target.files[0];
      if (!file) return;
      const reader = new FileReader();
      reader.onload = (evt) => {
          const text = evt.target.result;
          const proteins = parseFasta(text);
          setBatchStatus("processing");
          const results = processBatchList(proteins);
          setBatchResults(results);
          setProcessedCount(proteins.length);
          setBatchProgress(100);
          setBatchStatus("completed");
          downloadBatchCSV(results);
      };
      reader.readAsText(file);
  };

  // --- CHART DATA MEMO ---
  const chartData = useMemo(() => {
    if (!data || !data.length) return [];
    
    const ssDensity = calculateDensity(data.ssBonds, data.length, ssWindowSize, ssSmoothing);
    const nDensity = calculateDensity(data.nLinked, data.length, glycoWindowSize, glycoSmoothing);
    const oDensity = calculateDensity(data.oLinked, data.length, glycoWindowSize, glycoSmoothing);
    // Phospho uses Glyco params
    const pDensity = calculateDensity(data.phosphorylation, data.length, glycoWindowSize, glycoSmoothing);

    const step = Math.max(1, Math.ceil(data.length / 600)); 
    const points = [];

    for (let i = 0; i < data.length; i += step) {
      const isHighSS = ssDensity[i] >= ssThreshold;
      const isHighN = nDensity[i] >= glycoThreshold;
      const isOverlap = isHighSS && isHighN;

      points.push({
        pos: i,
        ss: ssDensity[i],
        nLinked: nDensity[i],
        oLinked: oDensity[i],
        phos: pDensity[i],
        overlapHeight: isOverlap ? Math.max(ssDensity[i], nDensity[i]) * 1.1 : 0 
      });
    }
    
    if (points.length > 0 && points[points.length - 1].pos < data.length) {
        points.push({ pos: data.length, ss: 0, nLinked: 0, oLinked: 0, phos: 0, overlapHeight: 0 });
    }
    
    return points;
  }, [data, ssWindowSize, ssThreshold, ssSmoothing, glycoWindowSize, glycoThreshold, glycoSmoothing]);

  const overlapRegions = useMemo(() => {
    if (!chartData.length) return [];
    const regions = [];
    let currentStart = null;

    chartData.forEach((pt) => {
      if (pt.overlapHeight > 0 && currentStart === null) {
        currentStart = pt.pos;
      } else if (pt.overlapHeight === 0 && currentStart !== null) {
        regions.push({ start: currentStart, end: pt.pos });
        currentStart = null;
      }
    });
    if (currentStart !== null) {
        regions.push({ start: currentStart, end: chartData[chartData.length - 1].pos });
    }
    return regions;
  }, [chartData]);

  const commonMargins = { top: 0, right: 0, bottom: 0, left: 0 };
  const xAxisHeight = 20; 

  return (
    <div className="min-h-screen bg-slate-50 text-slate-900 font-sans pb-20">
      <header className="bg-white border-b border-slate-200 sticky top-0 z-30">
        <div className="max-w-7xl mx-auto px-4 h-16 flex items-center justify-between">
          <div className="flex items-center gap-2">
            <div className="w-8 h-8 bg-indigo-600 rounded-lg flex items-center justify-center text-white font-bold shadow-lg shadow-indigo-200">
              P
            </div>
            <h1 className="text-xl font-bold tracking-tight text-slate-800">Protein<span className="text-indigo-600">CoLoc</span></h1>
          </div>
          
          <nav className="flex gap-1 bg-slate-100 p-1 rounded-lg">
            <button onClick={() => setView('calibration')} className={`px-4 py-1.5 rounded-md text-sm font-medium transition-all ${view === 'calibration' ? 'bg-white text-indigo-700 shadow-sm' : 'text-slate-500 hover:text-slate-700'}`}>1. Calibration</button>
            <button onClick={() => setView('analysis')} className={`px-4 py-1.5 rounded-md text-sm font-medium transition-all ${view === 'analysis' ? 'bg-white text-indigo-700 shadow-sm' : 'text-slate-500 hover:text-slate-700'}`}>2. Analysis</button>
            <button onClick={() => setView('database')} className={`px-4 py-1.5 rounded-md text-sm font-medium transition-all ${view === 'database' ? 'bg-white text-indigo-700 shadow-sm' : 'text-slate-500 hover:text-slate-700'}`}>3. Full Proteome</button>
          </nav>
        </div>
      </header>

      <main className="max-w-7xl mx-auto px-4 py-8 space-y-6">
        
        {view !== 'database' && (
            <section className="flex flex-col md:flex-row gap-4 items-end">
            <div className="flex-1 w-full">
                <label className="block text-xs font-semibold text-slate-500 uppercase tracking-wider mb-1">Target Gene Name</label>
                <div className="relative">
                <input 
                    type="text" 
                    value={geneName}
                    onChange={(e) => setGeneName(e.target.value)}
                    onKeyDown={(e) => e.key === 'Enter' && handleSearch()}
                    className="w-full pl-10 pr-4 py-3 border border-slate-300 rounded-xl focus:ring-2 focus:ring-indigo-500 focus:border-indigo-500 outline-none transition-shadow"
                    placeholder="e.g. FN1, ALB, INS..."
                />
                <Search className="w-5 h-5 text-slate-400 absolute left-3 top-3.5" />
                </div>
            </div>
            <button 
                onClick={handleSearch}
                disabled={loading}
                className="bg-slate-900 hover:bg-slate-800 text-white px-6 py-3 rounded-xl font-medium transition-colors flex items-center gap-2 disabled:opacity-50 shadow-lg shadow-slate-200"
            >
                {loading ? <Activity className="animate-spin w-5 h-5" /> : <Search className="w-5 h-5" />}
                Analyze Protein
            </button>
            </section>
        )}

        {error && (
          <div className="bg-red-50 border border-red-200 text-red-700 px-4 py-3 rounded-lg flex items-center gap-2">
            <AlertCircle className="w-5 h-5" />
            {error}
          </div>
        )}

        {view !== 'database' && data && chartData.length > 0 && (
          <>
            <div className="flex flex-wrap gap-4 items-center text-sm text-slate-600 bg-white p-4 rounded-xl border border-slate-200 shadow-sm">
              <span className="font-semibold text-slate-900">ID: {data.id}</span>
              <div className="w-px h-4 bg-slate-300 mx-2"></div>
              <span className="font-semibold text-slate-900">Species: {data.scientificName}</span>
              <div className="w-px h-4 bg-slate-300 mx-2"></div>
              <span className="font-semibold text-slate-900">Length: {data.length} aa</span>
              <div className="w-px h-4 bg-slate-300 mx-2"></div>
              <Badge color="blue" text={`${data.ssBonds.length} Disulfides`} />
              <Badge color="red" text={`${data.nLinked.length} N-Linked`} />
              <Badge color="purple" text={`${data.phosphorylation.length} Phospho`} />
            </div>

            {view === 'calibration' && (
              <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
                 <div className="lg:col-span-2 space-y-0">
                  <Card className="p-6 pb-6">
                    <div className="flex justify-between items-center mb-4">
                      <h2 className="text-lg font-bold text-slate-800 flex items-center gap-2">
                        <Settings className="w-5 h-5 text-indigo-500" />
                        Calibration Canvas
                      </h2>
                      <div className="text-xs text-slate-400">Sync ID: proteinView</div>
                    </div>
                    {/* Calibration Charts */}
                    <div className="flex items-end mb-1">
                        <div className="w-20 text-xs font-bold text-slate-400 text-right pr-3 pb-6">Density</div>
                        <div className="flex-1 h-64 relative border-b border-slate-200">
                            <ResponsiveContainer width="100%" height="100%">
                                <ComposedChart data={chartData} syncId="proteinView" margin={commonMargins}>
                                    <CartesianGrid strokeDasharray="3 3" vertical={false} stroke="#f1f5f9" />
                                    <XAxis dataKey="pos" type="number" domain={[0, data.length]} height={xAxisHeight} tick={{fontSize: 10, fill: '#94a3b8'}} axisLine={false} tickLine={false} />
                                    <YAxis hide domain={[0, 'auto']} />
                                    <Tooltip contentStyle={{ borderRadius: '8px', border: 'none', boxShadow: '0 4px 6px -1px rgb(0 0 0 / 0.1)' }} labelStyle={{ color: '#64748b', fontSize: '12px' }} />
                                    <Area type="monotone" dataKey="ss" stroke="#38bdf8" fill="#38bdf8" fillOpacity={0.2} strokeWidth={2} name="SS Density" />
                                    <Line type="monotone" dataKey="nLinked" stroke="#f43f5e" strokeWidth={2} dot={false} name="N-Linked Density" />
                                    <Line type="monotone" dataKey="oLinked" stroke="#eab308" strokeWidth={2} strokeDasharray="5 5" dot={false} name="O-Linked Density" />
                                    <Line type="monotone" dataKey="phos" stroke="#a855f7" strokeWidth={2} dot={false} name="Phospho Density" />
                                    <ReferenceLine y={ssThreshold} stroke="#38bdf8" strokeDasharray="3 3" opacity={0.5} />
                                    <ReferenceLine y={glycoThreshold} stroke="#f43f5e" strokeDasharray="3 3" opacity={0.5} />
                                </ComposedChart>
                            </ResponsiveContainer>
                        </div>
                    </div>
                    <div className="flex items-center mb-1">
                        <div className="w-20 text-xs font-bold text-sky-500 text-right pr-3">Raw SS</div>
                        <div className="flex-1 h-6 bg-slate-100 rounded overflow-hidden relative">
                            <ResponsiveContainer width="100%" height="100%">
                                <ComposedChart data={chartData} syncId="proteinView" margin={commonMargins}>
                                    <XAxis dataKey="pos" type="number" domain={[0, data.length]} height={xAxisHeight} hide />
                                    <YAxis hide />
                                    <Tooltip cursor={{ stroke: 'black', strokeWidth: 1 }} content={<></>} />
                                    {data.ssBonds.map((pos, i) => (
                                        <ReferenceLine key={`ss-${i}`} x={pos} stroke="#0ea5e9" strokeWidth={1} />
                                    ))}
                                </ComposedChart>
                            </ResponsiveContainer>
                        </div>
                    </div>
                    <div className="flex items-center mb-1">
                        <div className="w-20 text-xs font-bold text-rose-500 text-right pr-3">Raw N-Gly</div>
                        <div className="flex-1 h-6 bg-slate-100 rounded overflow-hidden relative">
                            <ResponsiveContainer width="100%" height="100%">
                                <ComposedChart data={chartData} syncId="proteinView" margin={commonMargins}>
                                    <XAxis dataKey="pos" type="number" domain={[0, data.length]} height={xAxisHeight} hide />
                                    <YAxis hide />
                                    <Tooltip cursor={{ stroke: 'black', strokeWidth: 1 }} content={<></>} />
                                    {data.nLinked.map((pos, i) => (
                                        <ReferenceLine key={`n-${i}`} x={pos} stroke="#f43f5e" strokeWidth={1} />
                                    ))}
                                </ComposedChart>
                            </ResponsiveContainer>
                        </div>
                    </div>
                     <div className="flex items-center mb-1">
                        <div className="w-20 text-xs font-bold text-yellow-600 text-right pr-3">Raw O-Gly</div>
                        <div className="flex-1 h-6 bg-slate-100 rounded overflow-hidden relative">
                            <ResponsiveContainer width="100%" height="100%">
                                <ComposedChart data={chartData} syncId="proteinView" margin={commonMargins}>
                                    <XAxis dataKey="pos" type="number" domain={[0, data.length]} height={xAxisHeight} hide />
                                    <YAxis hide />
                                    <Tooltip cursor={{ stroke: 'black', strokeWidth: 1 }} content={<></>} />
                                    {data.oLinked.map((pos, i) => (
                                        <ReferenceLine key={`o-${i}`} x={pos} stroke="#ca8a04" strokeWidth={1} />
                                    ))}
                                </ComposedChart>
                            </ResponsiveContainer>
                        </div>
                    </div>
                    <div className="flex items-center">
                        <div className="w-20 text-xs font-bold text-purple-600 text-right pr-3">Raw Phospho</div>
                        <div className="flex-1 h-6 bg-slate-100 rounded overflow-hidden relative">
                            <ResponsiveContainer width="100%" height="100%">
                                <ComposedChart data={chartData} syncId="proteinView" margin={commonMargins}>
                                    <XAxis dataKey="pos" type="number" domain={[0, data.length]} height={xAxisHeight} hide />
                                    <YAxis hide />
                                    <Tooltip cursor={{ stroke: 'black', strokeWidth: 1 }} content={<></>} />
                                    {data.phosphorylation.map((pos, i) => (
                                        <ReferenceLine key={`p-${i}`} x={pos} stroke="#9333ea" strokeWidth={1} />
                                    ))}
                                </ComposedChart>
                            </ResponsiveContainer>
                        </div>
                    </div>
                  </Card>
                </div>

                <div className="space-y-6">
                    <Card className="p-6 bg-white h-full overflow-y-auto max-h-[800px]">
                        <h3 className="font-bold text-slate-800 mb-4 flex items-center gap-2">
                             <Sliders className="w-4 h-4" /> Calibration Controls
                        </h3>
                        
                        <div className="mb-6 pb-6 border-b border-slate-100">
                            <h4 className="text-xs font-bold text-sky-600 uppercase tracking-wider mb-3 flex items-center gap-1">
                                <div className="w-2 h-2 bg-sky-500 rounded-full"></div> Disulfide (Cysteines)
                            </h4>
                            <ParameterControl label="Window Size" value={ssWindowSize} onChange={setSsWindowSize} min={10} max={200} step={5} unit="aa" colorClass="accent-sky-500" />
                            <ParameterControl label="Threshold" value={ssThreshold} onChange={setSsThreshold} min={1} max={10} step={0.5} unit="sites" colorClass="accent-sky-500" />
                            <ParameterControl label="Smoothing" value={ssSmoothing} onChange={setSsSmoothing} min={0} max={20} step={1} unit="px" colorClass="accent-sky-500" />
                        </div>

                        <div className="mb-6">
                            <h4 className="text-xs font-bold text-rose-600 uppercase tracking-wider mb-3 flex items-center gap-1">
                                <div className="w-2 h-2 bg-rose-500 rounded-full"></div> Glyco & Phospho
                            </h4>
                            <ParameterControl label="Window Size" value={glycoWindowSize} onChange={setGlycoWindowSize} min={10} max={200} step={5} unit="aa" colorClass="accent-rose-500" />
                            <ParameterControl label="Threshold" value={glycoThreshold} onChange={setGlycoThreshold} min={1} max={10} step={0.5} unit="sites" colorClass="accent-rose-500" />
                            <ParameterControl label="Smoothing" value={glycoSmoothing} onChange={setGlycoSmoothing} min={0} max={20} step={1} unit="px" colorClass="accent-rose-500" />
                        </div>

                        <div className="flex gap-2">
                            <button onClick={exportParameters} className="flex-1 bg-white border border-slate-300 hover:bg-slate-50 text-slate-700 py-3 rounded-lg font-medium transition-colors flex justify-center items-center gap-2 text-sm">
                                <Save className="w-4 h-4" /> Save Params
                            </button>
                            <button onClick={() => setView('analysis')} className="flex-1 bg-slate-900 hover:bg-slate-800 text-white py-3 rounded-lg font-medium transition-colors flex justify-center items-center gap-2">
                                Confirm <ChevronRight className="w-4 h-4" />
                            </button>
                        </div>
                    </Card>
                </div>
              </div>
            )}

            {view === 'analysis' && (
                <div className="space-y-6">
                    <Card className="p-6">
                        <div className="flex justify-between items-center mb-4">
                             <h3 className="font-bold text-slate-800">Analysis View</h3>
                             <label className="flex items-center gap-2 text-sm text-indigo-600 font-medium cursor-pointer hover:text-indigo-800">
                                <Upload className="w-4 h-4" /> Load Params (.txt)
                                <input type="file" className="hidden" accept=".txt" onChange={importParameters} />
                             </label>
                        </div>
                        <div className="h-96 w-full">
                            <ResponsiveContainer width="100%" height="100%">
                                <ComposedChart data={chartData} margin={{ top: 10, right: 30, left: 10, bottom: 30 }}>
                                    <CartesianGrid strokeDasharray="3 3" stroke="#f1f5f9" />
                                    <XAxis dataKey="pos" label={{ value: 'Residue Position', position: 'bottom', offset: 0 }} tick={{fontSize: 12}} />
                                    <YAxis label={{ value: 'Relative Density', angle: -90, position: 'insideLeft' }} />
                                    <Tooltip contentStyle={{ backgroundColor: 'rgba(255, 255, 255, 0.95)', borderRadius: '8px', boxShadow: '0 4px 12px rgba(0,0,0,0.1)', border: 'none' }} />
                                    {overlapRegions.map((region, i) => (
                                        <ReferenceArea key={i} x1={region.start} x2={region.end} fill="#86efac" fillOpacity={0.4} />
                                    ))}
                                    <Area type="monotone" dataKey="ss" stroke="#38bdf8" fill="#38bdf8" fillOpacity={0.1} strokeWidth={2} name="Disulfide" />
                                    <Line type="monotone" dataKey="nLinked" stroke="#f43f5e" strokeWidth={2} dot={false} name="N-Linked" />
                                    <Line type="monotone" dataKey="oLinked" stroke="#eab308" strokeWidth={2} strokeDasharray="5 5" dot={false} name="O-Linked" />
                                    <Line type="monotone" dataKey="phos" stroke="#a855f7" strokeWidth={2} dot={false} name="Phospho" />
                                </ComposedChart>
                            </ResponsiveContainer>
                        </div>
                    </Card>

                    <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                        <Card className="p-6">
                            <h3 className="font-bold text-slate-800 mb-4">Detected Regions of Interest</h3>
                            <div className="space-y-2 max-h-60 overflow-y-auto pr-2">
                                {overlapRegions.length === 0 ? (
                                    <p className="text-slate-500 text-sm italic">No significant co-localization zones found with current parameters.</p>
                                ) : (
                                    overlapRegions.map((reg, idx) => (
                                        <div key={idx} className="flex justify-between items-center p-3 bg-emerald-50 rounded-lg border border-emerald-100">
                                            <span className="font-mono text-emerald-800 font-medium">Region {idx + 1}</span>
                                            <span className="text-sm text-emerald-700">AA {reg.start} - {reg.end}</span>
                                            <span className="text-xs bg-white px-2 py-1 rounded text-emerald-600 border border-emerald-200 font-bold">
                                                {reg.end - reg.start} aa
                                            </span>
                                        </div>
                                    ))
                                )}
                            </div>
                        </Card>

                        <Card className="p-6 flex flex-col justify-center items-center text-center space-y-4">
                            <div className="w-12 h-12 bg-slate-100 rounded-full flex items-center justify-center text-slate-500">
                                <FileText className="w-6 h-6" />
                            </div>
                            <div>
                                <h3 className="font-bold text-slate-800">Export Results</h3>
                                <p className="text-sm text-slate-500 mt-1">Download a CSV report with densities and detected regions for {data.gene}.</p>
                            </div>
                            <button 
                                onClick={handleExport}
                                className="text-indigo-600 font-medium text-sm hover:underline flex items-center gap-2"
                            >
                                <Download className="w-4 h-4" /> Download CSV
                            </button>
                        </Card>
                    </div>
                </div>
            )}
          </>
        )}

        {view === 'database' && (
            <div className="space-y-6">
                <div className="bg-indigo-900 rounded-2xl p-8 text-white relative overflow-hidden">
                    <div className="relative z-10 max-w-3xl mx-auto">
                        <div className="flex justify-center items-center gap-4 mb-2">
                            <h2 className="text-2xl font-bold">Full Proteome Analysis</h2>
                            <label className="flex items-center gap-2 text-xs bg-indigo-700 hover:bg-indigo-600 px-3 py-1 rounded-full cursor-pointer transition-colors">
                                <Upload className="w-3 h-3" /> Load Params
                                <input type="file" className="hidden" accept=".txt" onChange={importParameters} />
                            </label>
                        </div>
                        <p className="text-indigo-200 mb-6 text-center">
                            This module will process thousands of proteins using your calibrated parameters: 
                            <br/>
                            <span className="font-mono text-xs bg-indigo-800 px-2 py-1 rounded ml-1">
                                SS[W:{ssWindowSize}/T:{ssThreshold}]
                            </span>
                            <span className="font-mono text-xs bg-indigo-800 px-2 py-1 rounded ml-1">
                                Glyco[W:{glycoWindowSize}/T:{glycoThreshold}]
                            </span>
                        </p>
                        
                        {batchStatus === 'idle' && (
                            <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                                {/* Card 1: Whole Organism */}
                                <div className="bg-indigo-800 bg-opacity-50 p-6 rounded-xl border border-indigo-700">
                                    <div className="flex items-center gap-3 mb-4">
                                        <Database className="w-6 h-6 text-indigo-300" />
                                        <h3 className="font-bold">UniProt (Reviewed)</h3>
                                    </div>
                                    <label className="block text-xs text-indigo-300 mb-1">Organism ID (TaxID)</label>
                                    <div className="flex gap-2">
                                        <input 
                                            type="text" 
                                            value={organismId} 
                                            onChange={(e) => setOrganismId(e.target.value)}
                                            className="bg-indigo-900 border border-indigo-600 rounded px-3 py-2 text-sm w-full focus:outline-none focus:border-white"
                                        />
                                        <button 
                                            onClick={runWholeOrganismAnalysis}
                                            className="bg-white text-indigo-900 px-4 py-2 rounded font-bold text-sm hover:bg-indigo-50 whitespace-nowrap flex items-center gap-2"
                                        >
                                            <Play className="w-3 h-3" /> Run
                                        </button>
                                    </div>
                                    <p className="text-xs text-indigo-400 mt-2">Default: 9606 (Homo sapiens)</p>
                                </div>

                                {/* Card 2: FASTA Upload */}
                                <div className="bg-indigo-800 bg-opacity-50 p-6 rounded-xl border border-indigo-700">
                                    <div className="flex items-center gap-3 mb-4">
                                        <FileText className="w-6 h-6 text-indigo-300" />
                                        <h3 className="font-bold">Custom FASTA</h3>
                                    </div>
                                    <label className="block text-xs text-indigo-300 mb-1">Upload .fasta file</label>
                                    <label className="flex items-center justify-center w-full h-10 px-4 transition bg-indigo-900 border border-indigo-600 border-dashed rounded appearance-none cursor-pointer hover:border-indigo-400 focus:outline-none">
                                        <span className="text-sm text-indigo-300">Select file...</span>
                                        <input type="file" accept=".fasta,.txt" className="hidden" onChange={handleFastaUpload} />
                                    </label>
                                    <p className="text-xs text-indigo-400 mt-2">Analysis based on motif search (N-X-S/T).</p>
                                </div>
                            </div>
                        )}

                        {batchStatus !== 'idle' && (
                            <div className="bg-white bg-opacity-10 p-6 rounded-xl text-center">
                                <div className="mb-4 flex justify-center">
                                    {batchStatus === 'processing' ? (
                                        <Loader className="w-12 h-12 animate-spin text-indigo-300" />
                                    ) : (
                                        <CheckCircle className="w-12 h-12 text-emerald-400" />
                                    )}
                                </div>
                                <h3 className="text-xl font-bold mb-2">
                                    {batchStatus === 'processing' ? 'Processing Proteome...' : 'Analysis Completed'}
                                </h3>
                                <p className="text-indigo-200 mb-4">Processed {processedCount} proteins</p>
                                
                                <div className="w-full bg-indigo-950 rounded-full h-2.5 mb-4 overflow-hidden">
                                    <div 
                                        className="bg-emerald-400 h-2.5 rounded-full transition-all duration-300" 
                                        style={{ width: `${batchProgress}%` }}
                                    ></div>
                                </div>
                                
                                {batchStatus === 'completed' && (
                                    <div className="mt-4">
                                        <button 
                                            onClick={() => setBatchStatus('idle')}
                                            className="text-sm underline text-indigo-300 hover:text-white"
                                        >
                                            Start New Analysis
                                        </button>
                                    </div>
                                )}
                            </div>
                        )}
                    </div>
                </div>
            </div>
        )}
      </main>
    </div>
  );
}
