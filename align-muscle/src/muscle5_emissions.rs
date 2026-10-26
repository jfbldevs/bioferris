// The BLOSUM joint-frequency record is kept at the source's published decimal
// precision, and exact zero/symmetry comparisons enforce table invariants.
#![allow(
    clippy::excessive_precision,
    clippy::float_cmp,
    clippy::unreadable_literal
)]

use core::any::TypeId;

use align_core::{Alphabet, Dna, Protein, Rna, Symbol};

const MAX_CANONICAL_SYMBOLS: usize = 20;

/// Canonical emission probabilities before optional HMM perturbation.
#[derive(Clone, Debug)]
pub(crate) struct EmissionModel {
    canonical_count: usize,
    joint: [[f64; MAX_CANONICAL_SYMBOLS]; MAX_CANONICAL_SYMBOLS],
    marginal: [f64; MAX_CANONICAL_SYMBOLS],
}

impl EmissionModel {
    pub(crate) fn for_alphabet<A: Muscle5Alphabet>() -> Self {
        base_emission_model::<A>()
    }

    pub(crate) const fn canonical_count(&self) -> usize {
        self.canonical_count
    }

    pub(crate) const fn joint_mut(
        &mut self,
    ) -> &mut [[f64; MAX_CANONICAL_SYMBOLS]; MAX_CANONICAL_SYMBOLS] {
        &mut self.joint
    }

    pub(crate) fn recompute_marginal(&mut self) {
        self.marginal.fill(0.0);
        for row in 0..self.canonical_count {
            self.marginal[row] = self.joint[row][..self.canonical_count].iter().sum();
        }
    }

    pub(crate) fn match_probability<A: Muscle5Alphabet>(
        &self,
        left: Symbol<A>,
        right: Symbol<A>,
    ) -> f64 {
        let left_mix = canonical_mixture(left, &self.marginal);
        let right_mix = canonical_mixture(right, &self.marginal);
        let mut probability = 0.0;
        for (first, first_weight) in left_mix[..self.canonical_count].iter().enumerate() {
            if *first_weight == 0.0 {
                continue;
            }
            for (second, second_weight) in right_mix[..self.canonical_count].iter().enumerate() {
                probability += first_weight * second_weight * self.joint[first][second];
            }
        }
        probability
    }

    pub(crate) fn insert_probability<A: Muscle5Alphabet>(&self, symbol: Symbol<A>) -> f64 {
        let mixture = canonical_mixture(symbol, &self.marginal);
        mixture[..self.canonical_count]
            .iter()
            .zip(&self.marginal[..self.canonical_count])
            .map(|(weight, probability)| weight * probability)
            .sum()
    }

    pub(crate) fn validate(&self) -> bool {
        let joint_sum = self.joint[..self.canonical_count]
            .iter()
            .flat_map(|row| row[..self.canonical_count].iter())
            .sum::<f64>();
        let marginal_sum = self.marginal[..self.canonical_count].iter().sum::<f64>();
        (joint_sum - 1.0).abs() <= 1.0e-12
            && (marginal_sum - 1.0).abs() <= 1.0e-12
            && self.joint[..self.canonical_count]
                .iter()
                .flat_map(|row| row[..self.canonical_count].iter())
                .all(|value| value.is_finite() && *value > 0.0)
            && self.marginal[..self.canonical_count]
                .iter()
                .all(|value| value.is_finite() && *value > 0.0)
    }
}

mod sealed {
    pub trait Sealed {}
}

impl sealed::Sealed for Protein {}
impl sealed::Sealed for Dna {}
impl sealed::Sealed for Rna {}

/// Typed alphabet supported by the MUSCLE 5 reconstruction.
///
/// This sealed marker is implemented for protein, DNA, and RNA.
pub trait Muscle5Alphabet: Alphabet + sealed::Sealed {}

impl Muscle5Alphabet for Protein {}

fn protein_emission_model() -> EmissionModel {
    let joint = BLOSUM62_JOINT;
    let mut marginal = [0.0; MAX_CANONICAL_SYMBOLS];
    for row in 0..20 {
        marginal[row] = joint[row].iter().sum();
    }
    EmissionModel {
        canonical_count: 20,
        joint,
        marginal,
    }
}

fn protein_mixture<A: Muscle5Alphabet>(
    symbol: Symbol<A>,
    background: &[f64; MAX_CANONICAL_SYMBOLS],
) -> [f64; MAX_CANONICAL_SYMBOLS] {
    let mut mixture = [0.0; MAX_CANONICAL_SYMBOLS];
    if symbol.alphabet_index() < 20 {
        mixture[symbol.alphabet_index()] = 1.0;
        return mixture;
    }
    match symbol.as_ascii() {
        b'B' => {
            mixture[2] = 0.5;
            mixture[3] = 0.5;
        }
        b'Z' => {
            mixture[5] = 0.5;
            mixture[6] = 0.5;
        }
        b'X' => mixture.copy_from_slice(background),
        b'J' => {
            mixture[9] = 0.5;
            mixture[10] = 0.5;
        }
        b'U' => mixture[4] = 1.0,
        b'O' => mixture[11] = 1.0,
        _ => unreachable!("validated protein symbol"),
    }
    mixture
}

impl Muscle5Alphabet for Dna {}
impl Muscle5Alphabet for Rna {}

pub(crate) fn canonical_word_index<A: Muscle5Alphabet>(symbol: Symbol<A>) -> Option<usize> {
    let canonical_count = if TypeId::of::<A>() == TypeId::of::<Protein>() {
        20
    } else {
        4
    };
    (symbol.alphabet_index() < canonical_count).then_some(symbol.alphabet_index())
}

pub(crate) fn super5_word_length<A: Muscle5Alphabet>() -> usize {
    if TypeId::of::<A>() == TypeId::of::<Protein>() {
        3
    } else {
        6
    }
}

fn base_emission_model<A: Muscle5Alphabet>() -> EmissionModel {
    if TypeId::of::<A>() == TypeId::of::<Protein>() {
        protein_emission_model()
    } else {
        nucleotide_emission_model()
    }
}

fn canonical_mixture<A: Muscle5Alphabet>(
    symbol: Symbol<A>,
    background: &[f64; MAX_CANONICAL_SYMBOLS],
) -> [f64; MAX_CANONICAL_SYMBOLS] {
    if TypeId::of::<A>() == TypeId::of::<Protein>() {
        protein_mixture(symbol, background)
    } else {
        nucleotide_mixture(symbol.as_ascii())
    }
}

fn nucleotide_emission_model() -> EmissionModel {
    let mut joint = [[0.0; MAX_CANONICAL_SYMBOLS]; MAX_CANONICAL_SYMBOLS];
    for (row, values) in joint[..4].iter_mut().enumerate() {
        for (column, value) in values[..4].iter_mut().enumerate() {
            *value = if row == column { 0.22 } else { 0.01 };
        }
    }
    let mut marginal = [0.0; MAX_CANONICAL_SYMBOLS];
    marginal[..4].fill(0.25);
    EmissionModel {
        canonical_count: 4,
        joint,
        marginal,
    }
}

fn nucleotide_mixture(ascii: u8) -> [f64; MAX_CANONICAL_SYMBOLS] {
    let mask: u8 = match ascii {
        b'A' => 0b0001,
        b'C' => 0b0010,
        b'G' => 0b0100,
        b'T' | b'U' => 0b1000,
        b'R' => 0b0101,
        b'Y' => 0b1010,
        b'S' => 0b0110,
        b'W' => 0b1001,
        b'K' => 0b1100,
        b'M' => 0b0011,
        b'B' => 0b1110,
        b'D' => 0b1101,
        b'H' => 0b1011,
        b'V' => 0b0111,
        b'N' => 0b1111,
        _ => unreachable!("validated nucleotide symbol"),
    };
    let count = f64::from(mask.count_ones());
    let mut mixture = [0.0; MAX_CANONICAL_SYMBOLS];
    for (index, weight) in mixture[..4].iter_mut().enumerate() {
        if mask & (1 << index) != 0 {
            *weight = 1.0 / count;
        }
    }
    mixture
}

/// Full-precision BLOSUM62 target frequencies from the NCBI public-domain record.
///
/// Canonical order is `ARNDCQEGHILKMFPSTWYV`.
pub(crate) const BLOSUM62_JOINT: [[f64; 20]; 20] = [
    [
        2.1497573378347484e-02,
        2.3470224274721213e-03,
        1.9493235258876179e-03,
        2.1674844853066858e-03,
        1.5903351423026848e-03,
        1.9242657898716525e-03,
        2.9879059292799641e-03,
        5.8158526388051033e-03,
        1.1076584657559144e-03,
        3.1880644746334580e-03,
        4.4186245468471547e-03,
        3.3466571942021082e-03,
        1.3412107617355408e-03,
        1.6360627863999076e-03,
        2.1568959784943114e-03,
        6.2524987419815400e-03,
        3.7180506975672363e-03,
        4.0281679108936688e-04,
        1.2999956675626666e-03,
        5.0679056444508912e-03,
    ],
    [
        2.3470224274721213e-03,
        1.7757465118386322e-02,
        1.9786027128591904e-03,
        1.5865480081162602e-03,
        3.9365984789376245e-04,
        2.4858611089731411e-03,
        2.6933867548771758e-03,
        1.7221140903704937e-03,
        1.2407382229440791e-03,
        1.2435878276496955e-03,
        2.4193952633248727e-03,
        6.2339060289407083e-03,
        8.0309461712520876e-04,
        9.3181986323789834e-04,
        9.5783034332718700e-04,
        2.2660898636037261e-03,
        1.7802796534180537e-03,
        2.6571979312581875e-04,
        9.2634607111251918e-04,
        1.5810185245264004e-03,
    ],
    [
        1.9493235258876179e-03,
        1.9786027128591904e-03,
        1.4140291972553610e-02,
        3.7201973506001745e-03,
        4.3845466068066216e-04,
        1.5304436972610567e-03,
        2.2097156829738759e-03,
        2.8591871815612977e-03,
        1.4301072616183181e-03,
        9.9437221166923172e-04,
        1.3690958423974782e-03,
        2.4402105140841090e-03,
        5.2943633069226512e-04,
        7.5004227978192801e-04,
        8.6016459857770028e-04,
        3.1466019144814608e-03,
        2.2360795375444384e-03,
        1.6159545671597605e-04,
        7.0048422794024819e-04,
        1.2014015528772706e-03,
    ],
    [
        2.1674844853066858e-03,
        1.5865480081162602e-03,
        3.7201973506001745e-03,
        2.1274574617480089e-02,
        3.9909227141697264e-04,
        1.6481246723433428e-03,
        4.9158017471929655e-03,
        2.5221102126636373e-03,
        9.5384849402143984e-04,
        1.2347404942429857e-03,
        1.5202051791453383e-03,
        2.4453087721980561e-03,
        4.6429229320514104e-04,
        7.6023722413111566e-04,
        1.2373315413524663e-03,
        2.8035127901697272e-03,
        1.8961512776990257e-03,
        1.6218020183662784e-04,
        5.9842263937853702e-04,
        1.3158365660538270e-03,
    ],
    [
        1.5903351423026848e-03,
        3.9365984789376245e-04,
        4.3845466068066216e-04,
        3.9909227141697264e-04,
        1.1931352277704348e-02,
        3.0937204045913537e-04,
        3.8338775043186374e-04,
        7.6951976030099293e-04,
        2.2976387481074697e-04,
        1.0956590131781735e-03,
        1.5682982157153873e-03,
        5.0124929379033781e-04,
        3.7717165634097634e-04,
        5.1389991547056834e-04,
        3.6111795849154795e-04,
        1.0432626586831986e-03,
        9.3041313726939057e-04,
        1.4474923964368156e-04,
        3.4603772624580643e-04,
        1.3606607271146112e-03,
    ],
    [
        1.9242657898716525e-03,
        2.4858611089731411e-03,
        1.5304436972610567e-03,
        1.6481246723433428e-03,
        3.0937204045913537e-04,
        7.3292255467189687e-03,
        3.5385780499965817e-03,
        1.3683038039160171e-03,
        1.0489026828741754e-03,
        8.9102936026571569e-04,
        1.6174411456311808e-03,
        3.0968229715707327e-03,
        7.3993258722701268e-04,
        5.4255147972143906e-04,
        8.4668181752066874e-04,
        1.8931125300036275e-03,
        1.3796838284921874e-03,
        2.2737931366728891e-04,
        6.7584155312457842e-04,
        1.1660966117775285e-03,
    ],
    [
        2.9879059292799641e-03,
        2.6933867548771758e-03,
        2.2097156829738759e-03,
        4.9158017471929655e-03,
        3.8338775043186374e-04,
        3.5385780499965817e-03,
        1.6133927472163669e-02,
        1.9380952488713059e-03,
        1.3667885452189439e-03,
        1.2192061706431622e-03,
        2.0030316026648431e-03,
        4.1322603720305197e-03,
        6.7909745467514783e-04,
        8.5179405867513139e-04,
        1.4216207127018586e-03,
        2.9539180653600089e-03,
        2.0493063257644955e-03,
        2.6488552587183780e-04,
        8.7044186256788659e-04,
        1.6987763526262680e-03,
    ],
    [
        5.8158526388051033e-03,
        1.7221140903704937e-03,
        2.8591871815612977e-03,
        2.5221102126636373e-03,
        7.6951976030099293e-04,
        1.3683038039160171e-03,
        1.9380952488713059e-03,
        3.7804346453413303e-02,
        9.5813607255887238e-04,
        1.3849118546156933e-03,
        2.0864716056392773e-03,
        2.5392537741810947e-03,
        7.3281559749652399e-04,
        1.1976708695723554e-03,
        1.3641171883713547e-03,
        3.8342830901664762e-03,
        2.1858459940987062e-03,
        4.0740829083805248e-04,
        8.3467413018106177e-04,
        1.8218235950233687e-03,
    ],
    [
        1.1076584657559144e-03,
        1.2407382229440791e-03,
        1.4301072616183181e-03,
        9.5384849402143984e-04,
        2.2976387481074697e-04,
        1.0489026828741754e-03,
        1.3667885452189439e-03,
        9.5813607255887238e-04,
        9.2802502369336622e-03,
        5.8089627083019206e-04,
        9.8696608463236094e-04,
        1.1873625842258938e-03,
        3.8264639620910225e-04,
        8.1041076335565583e-04,
        4.7770135861914477e-04,
        1.1052034635193162e-03,
        7.4371746073077327e-04,
        1.5168037757411286e-04,
        1.5213771111755425e-03,
        6.4882907765797669e-04,
    ],
    [
        3.1880644746334580e-03,
        1.2435878276496955e-03,
        9.9437221166923172e-04,
        1.2347404942429857e-03,
        1.0956590131781735e-03,
        8.9102936026571569e-04,
        1.2192061706431622e-03,
        1.3849118546156933e-03,
        5.8089627083019206e-04,
        1.8441526588740136e-02,
        1.1382470627796603e-02,
        1.5655862274689192e-03,
        2.5081290988482057e-03,
        3.0458868657559346e-03,
        1.0068164685944146e-03,
        1.7225081689171561e-03,
        2.6953622613315018e-03,
        3.6183761166072852e-04,
        1.3821121844492116e-03,
        1.1972663837662637e-02,
    ],
    [
        4.4186245468471547e-03,
        2.4193952633248727e-03,
        1.3690958423974782e-03,
        1.5202051791453383e-03,
        1.5682982157153873e-03,
        1.6174411456311808e-03,
        2.0030316026648431e-03,
        2.0864716056392773e-03,
        9.8696608463236094e-04,
        1.1382470627796603e-02,
        3.7141460156350926e-02,
        2.4634345023228079e-03,
        4.9293545515183088e-03,
        5.4151301166464015e-03,
        1.4146090399381900e-03,
        2.4277107072013821e-03,
        3.3238031308707055e-03,
        7.3206640617832933e-04,
        2.2096734692836624e-03,
        9.4786263030457313e-03,
    ],
    [
        3.3466571942021082e-03,
        6.2339060289407083e-03,
        2.4402105140841090e-03,
        2.4453087721980561e-03,
        5.0124929379033781e-04,
        3.0968229715707327e-03,
        4.1322603720305197e-03,
        2.5392537741810947e-03,
        1.1873625842258938e-03,
        1.5655862274689192e-03,
        2.4634345023228079e-03,
        1.6113385590544604e-02,
        9.0876633395557617e-04,
        9.4875149773685364e-04,
        1.5773020912564391e-03,
        3.1016069999481111e-03,
        2.3467014804084987e-03,
        2.7198500003555514e-04,
        9.9908866586876396e-04,
        1.9360424083099779e-03,
    ],
    [
        1.3412107617355408e-03,
        8.0309461712520876e-04,
        5.2943633069226512e-04,
        4.6429229320514104e-04,
        3.7717165634097634e-04,
        7.3993258722701268e-04,
        6.7909745467514783e-04,
        7.3281559749652399e-04,
        3.8264639620910225e-04,
        2.5081290988482057e-03,
        4.9293545515183088e-03,
        9.0876633395557617e-04,
        4.0477309321969848e-03,
        1.1901770463553603e-03,
        4.0824445213456919e-04,
        8.5603787638552766e-04,
        1.0095451907679563e-03,
        1.9872537223131380e-04,
        5.7145288352831449e-04,
        2.3123361470140736e-03,
    ],
    [
        1.6360627863999076e-03,
        9.3181986323789834e-04,
        7.5004227978192801e-04,
        7.6023722413111566e-04,
        5.1389991547056834e-04,
        5.4255147972143906e-04,
        8.5179405867513139e-04,
        1.1976708695723554e-03,
        8.1041076335565583e-04,
        3.0458868657559346e-03,
        5.4151301166464015e-03,
        9.4875149773685364e-04,
        1.1901770463553603e-03,
        1.8277684015431908e-02,
        5.2528021756783813e-04,
        1.1939618185901600e-03,
        1.1624184369750680e-03,
        8.4917468952377874e-04,
        4.2392005745634370e-03,
        2.5763052227920180e-03,
    ],
    [
        2.1568959784943114e-03,
        9.5783034332718700e-04,
        8.6016459857770028e-04,
        1.2373315413524663e-03,
        3.6111795849154795e-04,
        8.4668181752066874e-04,
        1.4216207127018586e-03,
        1.3641171883713547e-03,
        4.7770135861914477e-04,
        1.0068164685944146e-03,
        1.4146090399381900e-03,
        1.5773020912564391e-03,
        4.0824445213456919e-04,
        5.2528021756783813e-04,
        1.9066033679132538e-02,
        1.6662567934883051e-03,
        1.3511005665728870e-03,
        1.4152209821874487e-04,
        4.5224391125285910e-04,
        1.2451325046931832e-03,
    ],
    [
        6.2524987419815400e-03,
        2.2660898636037261e-03,
        3.1466019144814608e-03,
        2.8035127901697272e-03,
        1.0432626586831986e-03,
        1.8931125300036275e-03,
        2.9539180653600089e-03,
        3.8342830901664762e-03,
        1.1052034635193162e-03,
        1.7225081689171561e-03,
        2.4277107072013821e-03,
        3.1016069999481111e-03,
        8.5603787638552766e-04,
        1.1939618185901600e-03,
        1.6662567934883051e-03,
        1.2585947097159817e-02,
        4.7004857686835334e-03,
        2.8731729176487776e-04,
        1.0299846310599138e-03,
        2.3587292053265561e-03,
    ],
    [
        3.7180506975672363e-03,
        1.7802796534180537e-03,
        2.2360795375444384e-03,
        1.8961512776990257e-03,
        9.3041313726939057e-04,
        1.3796838284921874e-03,
        2.0493063257644955e-03,
        2.1858459940987062e-03,
        7.4371746073077327e-04,
        2.6953622613315018e-03,
        3.3238031308707055e-03,
        2.3467014804084987e-03,
        1.0095451907679563e-03,
        1.1624184369750680e-03,
        1.3511005665728870e-03,
        4.7004857686835334e-03,
        1.2514818886617953e-02,
        2.8575770858467209e-04,
        9.4161039895612720e-04,
        3.6402328079338207e-03,
    ],
    [
        4.0281679108936688e-04,
        2.6571979312581875e-04,
        1.6159545671597605e-04,
        1.6218020183662784e-04,
        1.4474923964368156e-04,
        2.2737931366728891e-04,
        2.6488552587183780e-04,
        4.0740829083805248e-04,
        1.5168037757411286e-04,
        3.6183761166072852e-04,
        7.3206640617832933e-04,
        2.7198500003555514e-04,
        1.9872537223131380e-04,
        8.4917468952377874e-04,
        1.4152209821874487e-04,
        2.8731729176487776e-04,
        2.8575770858467209e-04,
        6.4699301717154852e-03,
        8.8744160259272527e-04,
        3.5578318710317554e-04,
    ],
    [
        1.2999956675626666e-03,
        9.2634607111251918e-04,
        7.0048422794024819e-04,
        5.9842263937853702e-04,
        3.4603772624580643e-04,
        6.7584155312457842e-04,
        8.7044186256788659e-04,
        8.3467413018106177e-04,
        1.5213771111755425e-03,
        1.3821121844492116e-03,
        2.2096734692836624e-03,
        9.9908866586876396e-04,
        5.7145288352831449e-04,
        4.2392005745634370e-03,
        4.5224391125285910e-04,
        1.0299846310599138e-03,
        9.4161039895612720e-04,
        8.8744160259272527e-04,
        1.0246100213822419e-02,
        1.5489827890922993e-03,
    ],
    [
        5.0679056444508912e-03,
        1.5810185245264004e-03,
        1.2014015528772706e-03,
        1.3158365660538270e-03,
        1.3606607271146112e-03,
        1.1660966117775285e-03,
        1.6987763526262680e-03,
        1.8218235950233687e-03,
        6.4882907765797669e-04,
        1.1972663837662637e-02,
        9.4786263030457313e-03,
        1.9360424083099779e-03,
        2.3123361470140736e-03,
        2.5763052227920180e-03,
        1.2451325046931832e-03,
        2.3587292053265561e-03,
        3.6402328079338207e-03,
        3.5578318710317554e-04,
        1.5489827890922993e-03,
        1.9631915140537640e-02,
    ],
];

#[cfg(test)]
mod tests {
    #![allow(clippy::many_single_char_names)]

    use align_core::{Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence};

    use super::{EmissionModel, canonical_word_index};

    #[test]
    fn base_models_are_normalized_and_symmetric() {
        for model in [
            EmissionModel::for_alphabet::<Protein>(),
            EmissionModel::for_alphabet::<Dna>(),
            EmissionModel::for_alphabet::<Rna>(),
        ] {
            assert!(model.validate());
            for row in 0..model.canonical_count() {
                for column in 0..model.canonical_count() {
                    assert_eq!(model.joint[row][column], model.joint[column][row]);
                }
            }
        }
    }

    #[test]
    fn blosum_marginals_match_published_ncbi_values() {
        let model = EmissionModel::for_alphabet::<Protein>();
        let expected = [
            7.4216205067993410e-02,
            5.1614486141284638e-02,
            4.4645808512757915e-02,
            5.3626000838554413e-02,
            2.4687457167944848e-02,
            3.4259650591416023e-02,
            5.4311925684587502e-02,
            7.4146941452644999e-02,
            2.6212984805266227e-02,
            6.7917367618953756e-02,
            9.8907868497150955e-02,
            5.8155682303079680e-02,
            2.4990197579643110e-02,
            4.7418459742284751e-02,
            3.8538003320306206e-02,
            5.7229029476494421e-02,
            5.0891364550287033e-02,
            1.3029956129972148e-02,
            3.2281512313758580e-02,
            7.2919098205619245e-02,
        ];
        for (observed, expected) in model.marginal.iter().zip(expected) {
            assert!((observed - expected).abs() < 1.0e-15);
        }
    }

    #[test]
    fn ambiguity_emissions_are_explicit_mixtures() {
        let protein = EmissionModel::for_alphabet::<Protein>();
        let symbols = ProteinSequence::from_ascii("BDNXUO").unwrap();
        let b = symbols.get(0).unwrap();
        let d = symbols.get(1).unwrap();
        let x = symbols.get(3).unwrap();
        assert_eq!(
            protein.match_probability(b, d),
            0.5 * (protein.joint[2][3] + protein.joint[3][3])
        );
        assert!(protein.insert_probability(x) > 0.0);

        let dna = EmissionModel::for_alphabet::<Dna>();
        let n = DnaSequence::from_ascii("N").unwrap().get(0).unwrap();
        let a = DnaSequence::from_ascii("A").unwrap().get(0).unwrap();
        assert!((dna.match_probability(n, a) - 0.0625).abs() < 1.0e-15);

        let rna = EmissionModel::for_alphabet::<Rna>();
        let u = RnaSequence::from_ascii("U").unwrap().get(0).unwrap();
        assert!((rna.insert_probability(u) - 0.25).abs() < 1.0e-15);
    }

    #[test]
    fn word_index_rejects_ambiguities() {
        let canonical = ProteinSequence::from_ascii("A").unwrap().get(0).unwrap();
        let ambiguous = ProteinSequence::from_ascii("X").unwrap().get(0).unwrap();
        assert_eq!(canonical_word_index(canonical), Some(0));
        assert_eq!(canonical_word_index(ambiguous), None);
    }
}
