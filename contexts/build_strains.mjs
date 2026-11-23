// Laboratory strains of the microbial models as "contexts": the metabolic genes deleted or inactivated in their
// published genotype are switched off (no expression data and no task rescue: a deletion is a deletion).
//
//   node tools/contexts/build_strains.mjs        -> data/yeast-gem/contexts.json, data/ecoli-gem/contexts.json
//
// Genotypes: S. cerevisiae from SGD and Brachmann et al., Yeast 14:115 (1998) (BY series), Thomas & Rothstein,
// Cell 56:619 (1989) (W303), Entian & Kötter, Methods Microbiol 36:629 (2007) (CEN.PK); E. coli from Datsenko &
// Wanner, PNAS 97:6640 (2000) and Baba et al., Mol Syst Biol 2:2006.0008 (2006) (BW25113), Blattner et al., Science
// 277:1453 (1997) (MG1655), Taylor et al., J Microbiol Methods 17:233 (1993) (DH5α). Only lesions in genes of the model
// are applied; the others (recA, hsdR, rph, can1 alleles of non-metabolic genes…) are listed in the description.

const R = new URL("../../", import.meta.url).pathname;
const fs = await import("node:fs");
const { blockedReactions } = await import(R + "js/core/contexts.js");

const STRAINS = {
  "yeast-gem": [
    { id: "CEN.PK113-7D", desc: "MATa, prototrophic", genes: [], pairs: ["yeastO2_209", "yeastO2_28", "yeastO2_10", "yeastO2_05", "yeastO2_0"],
      note: "CEN.PK lineage, the reference strain of yeast physiology; Jouhten et al. used its sister CEN.PK113-1A (MATα)" },
    { id: "CEN.PK113-5D", desc: "MATa ura3-52", genes: ["URA3"], note: "needs uracil: use SC medium" },
    { id: "CEN.PK2-1C", desc: "MATa ura3-52 trp1-289 leu2-3,112 his3Δ1 MAL2-8c SUC2", genes: ["URA3", "TRP1", "LEU2", "HIS3"], note: "needs uracil, tryptophan, leucine and histidine: use SC medium" },
    { id: "BY4741", desc: "MATa his3Δ1 leu2Δ0 met15Δ0 ura3Δ0", genes: ["HIS3", "LEU2", "MET17", "URA3"], note: "S288C background; parent of the yeast deletion collection; needs histidine, leucine, methionine and uracil: use SC medium" },
    { id: "BY4742", desc: "MATα his3Δ1 leu2Δ0 lys2Δ0 ura3Δ0", genes: ["HIS3", "LEU2", "LYS2", "URA3"], note: "S288C background; needs histidine, leucine, lysine and uracil: use SC medium" },
    { id: "BY4743", desc: "MATa/α his3Δ1/his3Δ1 leu2Δ0/leu2Δ0 LYS2/lys2Δ0 MET15/met15Δ0 ura3Δ0/ura3Δ0", genes: ["HIS3", "LEU2", "URA3"], note: "diploid; lys2 and met15 are complemented; needs histidine, leucine and uracil: use SC medium" },
    { id: "W303-1A", desc: "MATa leu2-3,112 trp1-1 can1-100 ura3-1 ade2-1 his3-11,15", genes: ["LEU2", "TRP1", "CAN1", "URA3", "ADE2", "HIS3"], note: "can1-100: no arginine permease Can1; needs leucine, tryptophan, uracil, adenine and histidine: use SC medium" },
  ],
  "ecoli-gem": [
    { id: "MG1655", desc: "F⁻ λ⁻ ilvG⁻ rfb-50 rph-1", genes: [], note: "the strain of the model (K-12 reference); ilvG and rph are not active in the model either" },
    { id: "W3110", desc: "F⁻ λ⁻ IN(rrnD-rrnE)1 rph-1", genes: [], note: "K-12; metabolically the same genes as MG1655" },
    { id: "BW25113", desc: "Δ(araD-araB)567 Δ(rhaD-rhaB)568 ΔlacZ4787(::rrnB-3) hsdR514 rph-1", genes: ["araB", "araA", "araD", "rhaD", "rhaA", "rhaB", "lacZ"],
      pairs: ["ecoliZampieri", "ecoliGlcAer", "ecoliGlcAnaer"], note: "parent of the Keio deletion collection; cannot use arabinose, rhamnose or lactose" },
    { id: "DH5α", desc: "F⁻ φ80lacZΔM15 Δ(lacZYA-argF)U169 recA1 endA1 hsdR17 phoA supE44 thi-1 gyrA96 relA1", genes: ["lacZ", "lacY", "argF", "phoA"],
      note: "cloning strain; thi-1 (thiamine requirement) is not applied, its gene is not identified with certainty; argI still makes ornithine carbamoyltransferase; lacA is not in the model" },
  ],
};

for (const [slug, list] of Object.entries(STRAINS)) {
  const model = JSON.parse(fs.readFileSync(R + `data/${slug}/model.json`));
  const nG = model.genes.length, pos = new Map(model.genes.map((g, i) => [g[1], i]));
  const bits = (idx) => { const b = new Uint8Array(Math.ceil(nG / 8)); for (const i of idx) b[i >> 3] |= 1 << (i & 7); return Buffer.from(b).toString("base64"); };
  const contexts = list.map((s) => {
    const idx = s.genes.map((g) => { if (!pos.has(g)) throw new Error(`${slug} ${s.id}: no gene ${g}`); return pos.get(g); });
    const nBlocked = blockedReactions(model, new Set(s.genes)).length;
    return { id: "strain:" + s.id, name: s.id, group: "strain", src: "genotype", info: { genotype: s.desc, note: s.note || "" },
      off: bits(idx), nOff0: idx.length, rescued: [], nBlocked, failed: [], pairs: s.pairs || [] };
  });
  const doc = { meta: { sources: ["genotype"], threshold: "published genotype", license: "Genotypes from the literature (see docs)." },
    tasks: "none (genotypes)", genes: nG, taskNames: [], contexts };
  fs.writeFileSync(R + `data/${slug}/contexts.json`, JSON.stringify(doc));
  console.log(slug, contexts.map((c) => `${c.name}: ${c.nOff0} genes, ${c.nBlocked} reactions`).join("; "));
}
