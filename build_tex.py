#!/usr/bin/env python3
"""Convert PAPER.md to a REVTeX 4.2 manuscript (PRE format) and build the PDF.

Markdown here carries physics notation in Unicode (chi, sigma^2_v, <..>, t-hat, ...) rather
than in LaTeX, so the conversion is: protect code spans -> protect known math tokens ->
map Unicode -> escape remaining specials -> restore. Adjacent math fragments are merged into
single $...$ runs so the output reads naturally.

    python build_tex.py && latexmk -pdf paper.tex
"""
import re, sys

SRC, OUT = "PAPER.md", "paper.tex"

# ----------------------------------------------------------------- unicode
SUP = {'⁰':'0','¹':'1','²':'2','³':'3','⁴':'4','⁵':'5','⁶':'6','⁷':'7','⁸':'8','⁹':'9','⁻':'-'}
SUB = {'₁':'1','₂':'2','₃':'3','₄':'4'}
UNI = {
 '§':('\\S{}',2), '±':('\\pm',1), 'µ':('\\mu',1), '·':('\\cdot',1), '½':('1/2',1),
 'É':("\\'E",2), '×':('\\times',1), 'é':("\\'e",2), 'ö':('\\"o',2), 'ü':('\\"u',2),
 'ȧ':('\\dot{a}',1), 'Δ':('\\Delta',1), 'Σ':('\\sum',1), 'Φ':('\\Phi',1), 'Ψ':('\\Psi',1),
 'α':('\\alpha',1),'β':('\\beta',1),'δ':('\\delta',1),'ε':('\\varepsilon',1),'η':('\\eta',1),
 'θ':('\\theta',1),'λ':('\\lambda',1),'μ':('\\mu',1),'ξ':('\\xi',1),'π':('\\pi',1),
 'ρ':('\\rho',1),'σ':('\\sigma',1),'τ':('\\tau',1),'χ':('\\chi',1),'ψ':('\\psi',1),'ω':('\\omega',1),
 'ḃ':('\\dot{b}',1),'ḡ':('\\bar{g}',1),'Ṡ':('\\dot{S}',1),'ẋ':('\\dot{x}',1),'ẏ':('\\dot{y}',1),
 '–':('--',2), '—':('---',2), '‖':('\\|',1), '′':("'",1),
 '<':('<',1), '>':('>',1),
 '←':('\\leftarrow',1),'→':('\\rightarrow',1),'⇒':('\\Rightarrow',1),
 '∇':('\\nabla',1),'∈':('\\in',1),'−':('-',1),'∘':('\\circ',1),'∝':('\\propto',1),
 '∫':('\\int',1),'≈':('\\approx',1),'≠':('\\neq',1),'≡':('\\equiv',1),'≤':('\\le',1),
 '≥':('\\ge',1),'≳':('\\gtrsim',1),'⟨':('\\langle',1),'⟩':('\\rangle',1),
 '~':('\\sim',1), '|':('|',1),
}
ACC = {'\u0302':'hat', '\u0307':'dot', '\u0304':'bar'}

# math tokens written in plain text in the markdown, longest first
MATHTOK = [
 (r'⟨ΔX²⟩', r'\langle \Delta X^2\rangle'), (r'⟨ΔR²⟩', r'\langle \Delta R^2\rangle'),
 (r'⟨ΔS\|S⟩', r'\langle \Delta S|S\rangle'),
 (r'σ²_v', r'\sigma^2_v'), (r'kT_eff', r'k_BT_{\rm eff}'), (r'T_eff', r'T_{\rm eff}'),
 (r'ΔQ_null', r'\Delta Q_{\rm null}'), (r'n_cl', r'n_{\rm cl}'),
 (r'X_rel', r'X_{\rm rel}'), (r'S_NW', r'S_{\rm NW}'), (r'ρ_eq', r'\rho_{\rm eq}'),
 (r'ρ_ss', r'\rho_{\rm ss}'), (r'A_surr', r'A_{\rm surr}'), (r'k_B', r'k_B'),
 (r'F\^col', r'F^{\rm col}'), (r'F\^EHD', r'F^{\rm EHD}'), (r'L_x', r'L_x'), (r'L_y', r'L_y'),
 (r'N\^1\.00 ± 0\.07', r'N^{1.00\pm0.07}'),
 (r'R_turnover', r'R_{\rm turnover}'),
 (r'Ṡ_med', r'\dot S_{\rm med}'),
 (r'Σ_\{i∈species\}', r'\sum_{i\in\,{\rm species}}'),
 (r'⟨f⟩_row', r'\langle f\rangle_{\rm row}'),
 (r'⟨ψ⟩_s', r'\langle\psi\rangle_s'),
 (r'kT\b', r'k_BT'), (r'S\*', r'S^*'),
]
GREEK = 'αβγδεζηθικλμνξοπρστυφχψωΑΒΓΔΕΖΗΘΛΞΠΡΣΦΨΩ'
def to_math(t):
    """Map a short expression to pure math-mode LaTeX (used inside ^{...} and _{...})."""
    out = ''
    for c in t:
        if c in UNI: out += UNI[c][0]
        elif c in SUP: out += SUP[c]
        elif c in SUB: out += SUB[c]
        else: out += c
    return out
# X^{...}, X^-1/2, X^3.01 -- letters and Greek, ASCII caret
GEN_CARET = re.compile(r'([A-Za-z)\]'+GREEK+r'])\^(\{[^{}]*\}|[+−-]?[0-9]+(?:\.[0-9]+)?(?:/[0-9]+)?|[A-Za-z])')
GEN_SUP = re.compile(r'\b([A-Za-z])\^([0-9]+\.[0-9]+|[0-9]+)')
GEN_SUB = re.compile(r'(?<![A-Za-z\\])([A-Za-z'+GREEK+r'])_([A-Za-z]{1,4}|[0-9]{1,2})\b')

def sup_sub(t):
    # \x01 ... \x01 marks a fragment that must end up in math mode
    t = re.sub(r'([A-Za-z])?(['+''.join(SUP)+r']+)',
               lambda m: '\x01' + (m.group(1) or '') + '^{'
                         + ''.join(SUP[c] for c in m.group(2)) + '}\x01', t)
    t = re.sub('['+''.join(SUB)+']+',
               lambda m: '\x01_{'+''.join(SUB[c] for c in m.group(0))+'}\x01', t)
    return t

def uni_frag(text):
    """-> list of (latex, ismath)"""
    out, buf = [], ''
    i = 0
    while i < len(text):
        c = text[i]
        nxt = text[i+1] if i+1 < len(text) else ''
        if nxt in ACC:                                   # base + combining mark
            if buf: out.append((buf, 0)); buf = ''
            base = UNI[c][0] if c in UNI else c
            out.append(('\\%s{%s}' % (ACC[nxt], base), 1)); i += 2; continue
        if c == '\x01':
            if buf: out.append((buf, 0)); buf = ''
            j = text.index('\x01', i+1)
            out.append((text[i+1:j], 1)); i = j+1; continue
        if c == '√':
            if buf: out.append((buf, 0)); buf = ''
            m = re.match(r'√\(([^)]*)\)|√([A-Za-z0-9]+)', text[i:])
            arg = (m.group(1) or m.group(2)) if m else ''
            out.append(('\\sqrt{%s}' % arg, 1)); i += m.end() if m else 1; continue
        if c in UNI:
            if buf: out.append((buf, 0)); buf = ''
            out.append(UNI[c]); i += 1; continue
        buf += c; i += 1
    if buf: out.append((buf, 0))
    return out

ESC = {'&':r'\&', '%':r'\%', '#':r'\#', '$':r'\$', '~':r'\textasciitilde{}'}
def esc_text(t):
    t = t.replace('\\', r'\textbackslash{}')
    for k, v in ESC.items(): t = t.replace(k, v)
    t = t.replace('_', r'\_').replace('^', r'\textasciicircum{}')
    t = t.replace('{', r'\{').replace('}', r'\}')
    return t

def render(text):
    """Full inline conversion of one paragraph of markdown."""
    store = []
    def stash(latex):
        store.append(latex); return '\x00%d\x00' % (len(store)-1)
    # 0. raw LaTeX math already present in the markdown
    text = re.sub(r'\$([^$\n]+)\$', lambda m: stash('$' + m.group(1) + '$'), text)
    # 1. code spans
    text = re.sub(r'`([^`]+)`', lambda m: stash(r'\texttt{%s}' % m.group(1)
                  .replace('\\', r'\textbackslash{}').replace('_', r'\_')
                  .replace('&', r'\&').replace('%', r'\%').replace('#', r'\#')
                  .replace('$', r'\$').replace('{', r'\{').replace('}', r'\}')
                  .replace('^', r'\textasciicircum{}').replace('~', r'\textasciitilde{}')), text)
    # 1b. markdown escapes  \*  \|  \_  ->  literal char
    text = re.sub(r'\\([*|_\\[\]])', r'\1', text)
    # 2. citations
    def cite(m):
        body = m.group(1)
        ids = []
        for part in body.split(','):
            part = part.strip()
            if re.fullmatch(r'\d+', part): ids.append(part)
            elif re.fullmatch(r'\d+\s*[–-]\s*\d+', part):
                a, b = re.split(r'[–-]', part); ids += [str(x) for x in range(int(a), int(b)+1)]
            else: return m.group(0)
        return stash(r'\cite{%s}' % ','.join('r%s' % i for i in ids))
    text = re.sub(r'\[([0-9][0-9,\s–-]*)\]', cite, text)
    # 2b. bold vectors written as **r**_ij / **F**
    text = re.sub(r'\*\*([A-Za-z])\*\*_([A-Za-z0-9]+)',
                  lambda m: stash(r'$\mathbf{%s}_{%s}$' % (m.group(1), m.group(2))), text)
    text = re.sub(r'\*\*([A-Za-z])\*\*(?![A-Za-z])',
                  lambda m: stash(r'$\mathbf{%s}$' % m.group(1)), text)
    # 3. named math tokens, then generic sub/superscripts
    for pat, rep in MATHTOK:
        text = re.sub(pat, lambda m, r=rep: stash('$%s$' % r), text)
    def caret(m):
        base, expo = m.group(1), m.group(2)
        expo = expo[1:-1] if expo.startswith('{') else expo
        return stash('$%s^{%s}$' % (to_math(base), to_math(expo).replace('−', '-')))
    # standalone variable letters -> italic math (not sentence-initial, not 'A'/'I')
    text = re.sub(r'(?<=[\s(])([NQTLDSCXRZ])(?=[\s,.;:)])',
                  lambda m: stash('$%s$' % m.group(1)), text)
    text = GEN_CARET.sub(caret, text)
    text = GEN_SUP.sub(lambda m: stash('$%s^{%s}$' % (m.group(1), m.group(2))), text)
    def sub_(m):
        base, idx = to_math(m.group(1)), m.group(2)
        if len(idx) > 1 and idx.isalpha(): idx = '\\rm ' + idx
        return stash('$%s_{%s}$' % (base, idx))
    text = GEN_SUB.sub(sub_, text)
    # 4. unicode -> fragments, escaping the text ones, merging the math ones
    text = sup_sub(text)
    frags = uni_frag(text)
    out, mathbuf = [], []
    for lat, ism in frags:
        if ism == 1: mathbuf.append(lat)
        else:
            if mathbuf: out.append('$' + ' '.join(mathbuf) + '$'); mathbuf = []
            out.append(lat if ism == 2 else esc_text(lat))
    if mathbuf: out.append('$' + ' '.join(mathbuf) + '$')
    text = ''.join(out)
    # 4b. straight quotes -> TeX quotes
    text = re.sub(r'"([^"]*)"', lambda m: "``" + m.group(1) + "''", text)
    # 5. emphasis (after escaping, so ** survives)
    text = re.sub(r'\*\*(.+?)\*\*', r'\\textbf{\1}', text, flags=re.S)
    text = re.sub(r'(?<![\w*])\*([^*\n]+?)\*(?![\w*])', r'\\emph{\1}', text)
    # 6. restore
    text = re.sub(r'\x00(\d+)\x00', lambda m: store[int(m.group(1))], text)
    return text

# ----------------------------------------------------------------- blocks
EQ = {
"dx_i/dt = ξ_i(t) + (1/s_i) [ Σ_j F^col_ij + Σ_j F^EHD_ij ]":
 r"\frac{d\mathbf{x}_i}{dt}=\boldsymbol{\xi}_i(t)+\frac{1}{s_i}\Big[\sum_j\mathbf{F}^{\rm col}_{ij}+\sum_j\mathbf{F}^{\rm EHD}_{ij}\Big]",
"F^EHD_{i←j} = −α · l_j⁴ / (r² + l_j²)^{5/2} · **r**_ij ,    r ≤ 1":
 r"\mathbf{F}^{\rm EHD}_{i\leftarrow j}=-\alpha\,\frac{l_j^{4}}{(r^{2}+l_j^{2})^{5/2}}\,\mathbf{r}_{ij},\qquad r\le 1",
"a_i ≡ g(r; l_j) ,   a_j ≡ g(r; l_i)":
 r"a_i\equiv g(r;l_j),\qquad a_j\equiv g(r;l_i)",
"c_i = (1−χ)·ḡ + χ·a_i ,   c_j = (1−χ)·ḡ + χ·a_j\nF_i = −c_i **r**_ij ,     F_j = +c_j **r**_ij":
 r"""\begin{aligned}
c_i&=(1-\chi)\bar g+\chi a_i, &\qquad c_j&=(1-\chi)\bar g+\chi a_j,\\
\mathbf{F}_i&=-c_i\,\mathbf{r}_{ij}, &\qquad \mathbf{F}_j&=+c_j\,\mathbf{r}_{ij}.
\end{aligned}""",
"Q̇ = Σ_i **F**_i ∘ **ẋ**_i ,":
 r"\dot Q=\sum_i \mathbf{F}_i\circ\dot{\mathbf{x}}_i,",
"θ_i(t+dt) = arg⟨ e^{iθ_j} ⟩_{j ∈ N_i} + ω dt + η·U(−π, π)":
 r"\theta_i(t+dt)=\arg\big\langle e^{i\theta_j}\big\rangle_{j\in\mathcal{N}_i}+\omega\,dt+\eta\,U(-\pi,\pi)",
"T · EPR = χ J₁ + χ₂ J₂,    J₁ = Σ f₁ · v,   J₂ = Σ f₂ · v,":
 r"T\,\dot S=\chi J_1+\chi_2 J_2,\qquad J_1=\sum \mathbf{f}_1\!\cdot\!\mathbf{v},\qquad J_2=\sum \mathbf{f}_2\!\cdot\!\mathbf{v},",
"χ ⟨F_a∘ẋ⟩ = χ Σ_i (1/s_i) ⟨F_a,i · F_i⟩ + χ Σ_i D_i ⟨∇_i · F_a,i⟩ ,":
 r"\chi\,\langle\mathbf{F}_a\circ\dot{\mathbf{x}}\rangle=\chi\sum_i\frac{1}{s_i}\langle\mathbf{F}_{a,i}\cdot\mathbf{F}_i\rangle+\chi\sum_i D_i\langle\nabla_i\cdot\mathbf{F}_{a,i}\rangle,",
"S_NW = ∫ (E − I + A·C) dt":
 r"S_{\rm NW}=\int (E-I+A\cdot C)\,dt",
}

CAPS = {
1:"The simulation campaign.",
2:"Structural observables against $\\chi$ at paper scale ($N=4000$, three seeds), with the monodisperse reference for comparison.",
3:"Excess dissipation per particle against $\\chi$ at $N=4000$ from the configurational estimator on re-equilibrated configurations (three seeds), with its ratio to $\\chi^2$ and the contact-screening fraction. Parenthesised uncertainties apply to the last digits.",
4:"Condensate shedding rate and fragment reabsorption probability against $\\chi$ from the dense continuations, and the fragment count their ratio predicts against the count observed.",
5:"Complete source--sink budget of the fragment population per snapshot (seed 1 at each $\\chi$), closing to within 0.26 events per snapshot; each channel pair balances separately.",
6:"Attribution of the monodisperse-to-bidisperse difference between polydispersity and reciprocity, at paper scale.",
7:"Box scaling at fixed density $6.944$ particles per unit area, $\\chi=1.5$, all points continued to a late-time plateau.",
8:"Replicate at the source specification ($N=22{,}000$, $22.7\\%$ type-I): largest-cluster fraction by seed and observation window.",
9:"Becker--D\\\"oring size drift per interval against cluster size at $\\chi=1.5$, from monomer and dimer exchange events only ($\\Delta t=50$, three seeds).",
10:"Committor of fragment clusters at $\\chi=1.5$ from 72 relaunches of $\\hat t=3000$ with fresh noise, absorption into the condensate counted as the condensed outcome.",
11:"Horizon dependence of the half-committor crossing and of the fraction of fates decided near the crossover.",
12:"Density dependence at fixed $N=4000$ and $\\chi=1.5$: the condensate is strongly density-dependent, the fragment scale much less so.",
13:"A working taxonomy of commonly invoked variational tiers and the experimental signatures usually associated with them.",
14:"Differential cross-response between two nonreciprocal channels of \\emph{identical} structure, at the operating point $(\\chi,\\chi_2)=(0.5,0.5)$.",
15:"The same measurement with a structurally \\emph{dissimilar} second channel. The pattern inverts, which is why the reciprocity reading is withdrawn.",
16:"Single-particle effective temperature against lag and $\\chi$ from random-sign forcing with a shared noise stream (twelve twins per $\\chi$).",
17:"For the record: displacement scaling and effective temperature on species coordinates at $\\chi=0$ and $\\chi=1.5$, raw and with the system centre-of-mass drift removed.",
18:"\\emph{Chlamydomonas} axonemes by ATP concentration: beat frequency, circulation $z$-scores, and enclosed area per cycle.",
19:"Circulation at a matched level of description: single tracked colloidal clusters against single axonemes.",
20:"Vicsek-type ensembles in four variants: polar order against circulation.",
21:"Sampling floor of the Hodge cyclic fraction on synthetic perfectly transitive groups, raw counts against log-odds.",
22:"Hodge decomposition of four published sociomatrices, each against a matched Bradley--Terry floor.",
23:"Pre-registered predictions against outcomes. Thresholds were committed in \\texttt{EXPERIMENT.md} before any $\\chi$ run.",
}

FIGS = {  # anchor heading -> (file, label, caption)
"### 2.2 The reciprocity mixing parameter χ": ("fig1_construction",
 "fig:construction",
 "\\textbf{The $\\chi$ construction and its validation.} (a) Force coefficients $c_i$ (red) and "
 "$c_j$ (blue) against pair separation at $\\chi=0,\\,0.5,\\,1$. At $\\chi=0$ the two coincide, "
 "restoring Newton's third law exactly. (b) The symmetric sector (blue, five superposed curves at "
 "$\\chi=0,\\,0.25,\\,0.5,\\,1,\\,1.5$) is invariant in $\\chi$, while the antisymmetric sector "
 "(red) is exactly linear in it. (c) The antisymmetric-to-symmetric ratio for four values of the "
 "coupling strength $\\hat\\alpha$; the curves are indistinguishable, showing that $\\hat\\alpha$ "
 "cancels identically from the ratio and therefore cannot serve as a reciprocity knob."),
"### 3.2 Nonreciprocity produces the fragmented morphology termed arrested coarsening": ("fig2_structure",
 "fig:structure",
 "\\textbf{Structure against $\\chi$ at paper scale ($N=4000$, three seeds).} (a) Cluster count, "
 "log axis, showing the non-monotonic dip at $\\chi=0.25$ to a value five times \\emph{below} the "
 "reciprocal case before rising by two orders of magnitude; the dashed line is the monodisperse "
 "reference. (b) Largest-cluster fraction. (c) Activity $\\sigma^2_v$ (orange, log axis) against "
 "Newman modularity $Q$ (blue, right axis) on a scale spanning only $0.86$--$0.95$: the currents "
 "rise by two orders of magnitude while $Q$ barely moves. Error bars are the standard error over seeds."),
"### 3.8 No characteristic cluster size is selected; the small-cluster state is transient": ("fig3_scaling",
 "fig:scaling",
 "\\textbf{Box scaling at fixed density.} (a) Largest cluster against $N$ for three converged "
 "system sizes, with the phase-separation ($N^1$) and finite-characteristic-size ($N^0$) "
 "expectations. The measured exponent is $1.00\\pm0.07$. (b) The largest-cluster fraction is flat "
 "across a fourfold range of $N$, which is the same statement without a fit."),
"### 3.9 The selected scale is an event-rate crossover, not a stable cluster size": ("fig4_rates",
 "fig:rates",
 "\\textbf{The event-rate crossover.} (a) Per-cluster split and merge rates against cluster size, "
 "measured at $\\Delta t=50$ from equilibrated configurations ($105{,}708$ cluster observations). "
 "The rates cross between the size-5 and size-10 bins. (b) Their difference, showing the sign "
 "change. Splitting dominates below the crossing, merging above it."),
"### 3.10 The limits of three commonly used nonequilibrium diagnostics": ("fig5_thermo",
 "fig:thermo",
 "\\textbf{Dissipation and response.} (a) Excess dissipation per particle against $\\chi$, each "
 "point probed from a configuration equilibrated at that $\\chi$, with the quadratic form. "
 "(b) Symmetric and antisymmetric parts of the differential cross-response between two "
 "nonreciprocal channels, against measurement window. For channels of identical structure the "
 "symmetric part is consistent with zero at every window while the antisymmetric part is stable; "
 "the structurally dissimilar control is plotted alongside on the same axes, where the pattern "
 "inverts. (c) Effective temperature against measurement window: flat at $\\chi=0$, where it "
 "calibrates to $k_BT$, and growing as $t^{0.89}$ at $\\chi=1.5$ with no plateau. Inset: "
 "mean-squared displacement of the system centre of mass, near-ballistic at $\\chi=1.5$ over these "
 "windows ($t^{1.96}$) and diffusive at $\\chi=0$ ($t^{1.01}$)."),
"### 3.11 A biological comparison: does irreversibility survive coarse-graining?": ("fig6_circulation",
 "fig:circulation",
 "\\textbf{Circulation at matched level of description.} (a) Median $|z|$ of the signed area rate "
 "for a single tracked colloidal cluster at three values of $\\chi$, against a single "
 "\\emph{Chlamydomonas} axoneme, using the identical estimator, shape-descriptor reduction to two "
 "modes, and phase-randomised surrogate protocol. (b) The fraction of objects exceeding $|z|=2$. "
 "No circulation is detectable in the tracked colloidal clusters; axonemes show it in 92\\% of cases."),
}

def split_row(line):
    """Split a markdown table row, honouring \\| as a literal pipe."""
    line = line.strip()
    if line.startswith('|'): line = line[1:]
    if line.endswith('|'): line = line[:-1]
    cells, cur, i = [], '', 0
    while i < len(line):
        if line[i] == '\\' and i+1 < len(line) and line[i+1] == '|':
            cur += '|'; i += 2
        elif line[i] == '|':
            cells.append(cur); cur = ''; i += 1
        else:
            cur += line[i]; i += 1
    cells.append(cur)
    return [c.strip() for c in cells]

def emit_table(hdr, align, rows, idx):
    ncol = len(hdr)
    spec = ''.join('r' if a.endswith(':') and not a.startswith(':') else
                   ('c' if a.startswith(':') and a.endswith(':') else 'l') for a in align)
    spec = spec[:ncol].ljust(ncol, 'l')
    # long prose cells overflow an l/c/r spec; switch to proportional p{} columns
    widths = [max([len(render(hdr[c]))] + [len(render(r[c])) for r in rows if c < len(r)])
              for c in range(ncol)]
    if sum(widths) > 95:
        # weight by sqrt of content length: linear weighting starves the short columns
        import math
        wts = [math.sqrt(max(w, 1)) for w in widths]
        tot = float(sum(wts))
        avail = 15.6 - 0.45 * ncol      # cm of text width less inter-column separation
        spec = ''.join('p{%.2fcm}' % max(1.0, avail * w / tot) for w in wts)
    L = ['\\begin{table}[htbp]', '\\caption{%s}' % CAPS.get(idx, ''),
         '\\label{tab:%d}' % idx, '\\begin{ruledtabular}',
         '\\begin{tabular}{%s}' % spec]
    ragged = 'p{' in spec
    if ragged: L.append('\\raggedright')
    L.append(' & '.join(render(h) for h in hdr) + r' \\')
    L.append('\\hline')
    for r in rows:
        r = (r + ['']*ncol)[:ncol]
        L.append(' & '.join(render(c) for c in r) + r' \\')
    L += ['\\end{tabular}', '\\end{ruledtabular}', '\\end{table}']
    return '\n'.join(L)

def emit_fig(key):
    f, lab, cap = FIGS[key]
    return ('\\begin{figure}[htbp]\n\\centering\n'
            '\\includegraphics[width=\\linewidth]{%s.png}\n'
            '\\caption{%s}\n\\label{%s}\n\\end{figure}' % (f, cap, lab))

def convert(md):
    lines = md.split('\n')
    body, i, tno = [], 0, 0
    in_refs = False
    refs = []
    while i < len(lines):
        L = lines[i]
        # ---- headings
        if L.startswith('## '):
            t = L[3:].strip()
            if t.startswith('References'): in_refs = True; i += 1; continue
            if t.startswith('Abstract') or t.startswith('Figures'):
                # skip abstract (handled separately) / drop the figure-caption list
                i += 1
                while i < len(lines) and not lines[i].startswith('## '): i += 1
                continue
            t = re.sub(r'^\d+\.\s*', '', t)
            body.append('\n\\section{%s}' % render(t)); i += 1; continue
        if L.startswith('### '):
            key = L.rstrip()
            t = re.sub(r'^###\s*\d+\.\d+\s*', '', L).strip()
            body.append('\n\\subsection{%s}' % render(t))
            if key in FIGS: body.append(emit_fig(key))
            i += 1; continue
        if in_refs:
            m = re.match(r'^\[(\d+)\]\s*(.*)', L)
            if m:
                num, txt = m.group(1), m.group(2)
                j = i+1
                while j < len(lines) and lines[j].strip() and not lines[j].startswith('['):
                    txt += ' ' + lines[j].strip(); j += 1
                refs.append((num, render(txt))); i = j; continue
            i += 1; continue
        # ---- tables
        if L.strip().startswith('|') and i+1 < len(lines) and re.match(r'^\s*\|[\s:|-]+\|\s*$', lines[i+1]):
            tno += 1
            hdr = split_row(L); align = split_row(lines[i+1]); rows = []
            j = i+2
            while j < len(lines) and lines[j].strip().startswith('|'):
                rows.append(split_row(lines[j])); j += 1
            body.append(emit_table(hdr, align, rows, tno)); i = j; continue
        # ---- equations (4-space indent, not a list continuation)
        if re.match(r'^    \S', L) and not body_is_list(body):
            blk = []
            j = i
            while j < len(lines) and re.match(r'^    \S', lines[j]):
                blk.append(lines[j][4:].rstrip()); j += 1
            key = '\n'.join(blk)
            if key in EQ:
                body.append('\\begin{equation}\n%s\n\\end{equation}' % EQ[key]); i = j; continue
            body.append('\\begin{quote}\\ttfamily\n' + '\\\\\n'.join(esc_text(b) for b in blk) + '\n\\end{quote}')
            i = j; continue
        # ---- lists
        if re.match(r'^\d+\.\s', L) or re.match(r'^-\s', L):
            ordered = bool(re.match(r'^\d+\.\s', L))
            items, j = [], i
            while j < len(lines):
                m = re.match(r'^(\d+\.|-)\s+(.*)', lines[j])
                if m:
                    txt = m.group(2); j += 1
                    while j < len(lines) and lines[j].strip() and not re.match(r'^(\d+\.|-)\s', lines[j]) \
                          and not lines[j].startswith('#') and not lines[j].strip().startswith('|'):
                        txt += ' ' + lines[j].strip(); j += 1
                    items.append(txt)
                elif lines[j].strip() == '' and j+1 < len(lines) and re.match(r'^(\d+\.|-)\s', lines[j+1]):
                    j += 1
                else: break
            env = 'enumerate' if ordered else 'itemize'
            body.append('\\begin{%s}' % env)
            for it in items: body.append('\\item %s' % render(it))
            body.append('\\end{%s}' % env)
            i = j; continue
        # ---- horizontal rule / blank
        if L.strip() in ('---', ''): body.append(''); i += 1; continue
        # ---- paragraph
        para = [L]; j = i+1
        while j < len(lines) and lines[j].strip() and not lines[j].startswith('#') \
              and not lines[j].strip().startswith('|') and not re.match(r'^(\d+\.|-)\s', lines[j]) \
              and not re.match(r'^    \S', lines[j]) and lines[j].strip() != '---':
            para.append(lines[j]); j += 1
        body.append(render(' '.join(x.strip() for x in para))); i = j
    return '\n'.join(body), refs

def body_is_list(body):
    for b in reversed(body[-3:]):
        if b.strip(): return b.strip().startswith('\\item') or b.strip() == '\\begin{enumerate}'
    return False

PREAMBLE = r"""\documentclass[aps,pre,preprint,superscriptaddress,nofootinbib,longbibliography]{revtex4-2}
\usepackage{amsmath,amssymb}
\usepackage{graphicx}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{bm}
\usepackage{xcolor}
\usepackage[colorlinks=true,linkcolor=blue,citecolor=blue,urlcolor=blue]{hyperref}
\usepackage{textcomp}
\usepackage[htt]{hyphenat}
\setlength{\emergencystretch}{4em}
\tolerance=2000
\hyphenpenalty=1000
\begin{document}

\title{Isolating the antisymmetric sector of a nonreciprocal colloidal model:
kinetic unjamming, transient arrested coarsening, and irreversibility
without coarse-grained circulation}

\author{Martin G. Frasch}
\email{mfrasch@uw.edu}
\affiliation{Institute on Human Development and Disability, University of Washington,
Seattle, Washington 98195, USA}
\affiliation{Health Stream Analytics LLC, Seattle, Washington, USA}

\date{\today}

\begin{abstract}
%(ABSTRACT)s
\end{abstract}

\keywords{nonreciprocal interactions, active matter, arrested coarsening, entropy production,
kinetic arrest, coarse-graining}

\maketitle
"""

def main():
    md = open(SRC).read()
    # ---- pull the abstract out
    a0 = md.index('## Abstract') + len('## Abstract')
    a1 = md.index('**Keywords:**')
    abstract = md[a0:a1].strip()
    abstract = '\n\n'.join(render(' '.join(p.split()))
                           for p in re.split(r'\n\s*\n', abstract) if p.strip())
    # ---- body starts at the Introduction
    body_md = md[md.index('## 1. Introduction'):]
    body, refs = convert(body_md)
    bib = ['\n\\begin{thebibliography}{99}']
    for num, txt in refs:
        bib.append('\\bibitem{r%s} %s' % (num, txt))
    bib.append('\\end{thebibliography}')
    tex = (PREAMBLE % {'ABSTRACT': abstract}) + body + '\n'.join(bib) + '\n\n\\end{document}\n'
    tex = re.sub(r'\n{3,}', '\n\n', tex)
    open(OUT, 'w').write(tex)
    print('wrote %s  (%d lines, %d refs)' % (OUT, tex.count('\n'), len(refs)))

if __name__ == '__main__':
    main()
