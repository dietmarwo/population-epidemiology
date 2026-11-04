//! Exponential Ziggurat adapted from Julia Random; see THIRD_PARTY_NOTICES.md.
use super::{
    Bail, Draws, numeric,
    ziggurat_tables::{FE, KE, R, WE},
};
use rand_core::RngCore;

pub(super) fn exponential(rng: &mut impl RngCore, draws: &mut Draws) -> Result<f64, Bail> {
    exponential_from_words(|retry| {
        if retry {
            draws.retry()?;
        }
        draws.word(rng)
    })
}

/// Draw a unit exponential. The caller accounts for each word and can poll
/// for cancellation before a redraw following rejection (the flag is true).
#[inline]
pub fn exponential_from_words<E>(
    mut next_word: impl FnMut(bool) -> Result<u64, E>,
) -> Result<f64, E> {
    let mut retry = false;
    loop {
        let word = next_word(retry)? & 0x000f_ffff_ffff_ffff;
        // Exclude exact zero, as in the legacy inverse sampler.
        if word != 0 {
            let index = (word & 255) as usize;
            let x = word as f64 * WE[index];
            if word < KE[index] {
                return Ok(x);
            }
            let u = crate::sampling::uniform_closed_open(next_word(false)?);
            if index == 0 {
                return Ok(R - libm::log1p(-u));
            }
            let y = (FE[index - 1] - FE[index]) * u + FE[index];
            if y < libm::exp(-x) {
                return Ok(x);
            }
        }
        retry = true;
    }
}

pub(super) fn wait(total: f64, rng: &mut impl RngCore, draws: &mut Draws) -> Result<f64, Bail> {
    let value = exponential(rng, draws)? / total;
    if !value.is_finite() || value <= 0.0 {
        return Err(numeric("Ziggurat wait is outside the positive finite range").into());
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        sampling,
        simulation::{Control, Poller},
        ssa::Work,
    };
    use std::time::Instant;

    #[test]
    fn exponential_million_sample_cdf_gate() {
        let c = Control::default();
        let mut p = Poller::new(&c, Instant::now(), None);
        let mut w = Work::default();
        let mut d = Draws {
            limit: 10_000_000,
            work: &mut w,
            poller: &mut p,
            fast: false,
        };
        let mut rng = sampling::stream(4242, 778);
        let mut bins = [0u64; 91];
        let mut mean = 0.0;
        let mut second = 0.0;
        const N: u64 = 1_000_000;
        for _ in 0..N {
            let x = exponential(&mut rng, &mut d).ok().unwrap();
            assert!(x > 0.0 && x.is_finite());
            bins[(x * 10.0).floor().min(90.0) as usize] += 1;
            mean += x / N as f64;
            second += x * x / N as f64;
        }
        let bound = (libm::log(2.0 * 100.0 / 0.001) / (2.0 * N as f64)).sqrt();
        let mut cumulative = 0;
        for (k, count) in bins.iter().take(90).enumerate() {
            cumulative += count;
            let expected = -libm::expm1(-((k + 1) as f64) / 10.0);
            assert!((cumulative as f64 / N as f64 - expected).abs() <= bound);
        }
        eprintln!(
            "Ziggurat N={N} mean={mean} second_moment={second} rng_words={}",
            d.work.rng_words
        );
    }

    struct Words {
        words: Vec<u64>,
        position: usize,
    }
    impl RngCore for Words {
        fn next_u64(&mut self) -> u64 {
            let w = self.words[self.position % self.words.len()];
            self.position += 1;
            w
        }
        fn next_u32(&mut self) -> u32 {
            self.next_u64() as u32
        }
        fn fill_bytes(&mut self, _: &mut [u8]) {
            unreachable!()
        }
        fn try_fill_bytes(&mut self, _: &mut [u8]) -> Result<(), rand_core::Error> {
            unreachable!()
        }
    }
    #[test]
    fn zero_tail_and_rejection_branches_respect_budgets_and_stops() {
        let c = Control::default();
        let mut p = Poller::new(&c, Instant::now(), None);
        let mut w = Work::default();
        let mut d = Draws {
            limit: 4,
            work: &mut w,
            poller: &mut p,
            fast: false,
        };
        let mut rng = Words {
            words: vec![0],
            position: 0,
        };
        assert!(matches!(
            exponential(&mut rng, &mut d),
            Err(Bail::SamplingLimit)
        ));
        assert_eq!(d.work.rng_words, 4);
        d.work.rng_words = 0;
        rng.words = vec![0x000f_ffff_ffff_ff00, 1u64 << 63];
        rng.position = 0;
        let x = exponential(&mut rng, &mut d).ok().unwrap();
        assert_eq!(x, R + libm::log(2.0));
        assert_eq!(d.work.rng_words, 2);
        // Layer one always enters the slow branch; a high ordinate rejects.
        d.work.rng_words = 0;
        rng.words = vec![0x000f_ffff_ffff_ff01, u64::MAX];
        rng.position = 0;
        assert!(matches!(
            exponential(&mut rng, &mut d),
            Err(Bail::SamplingLimit)
        ));
        assert_eq!(d.work.rng_words, 4);
        c.cancel();
        assert!(matches!(
            exponential(&mut rng, &mut d),
            Err(Bail::Stopped(_))
        ));
    }
}
