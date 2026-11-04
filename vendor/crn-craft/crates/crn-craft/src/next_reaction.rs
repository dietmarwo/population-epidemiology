//! Indexed clock heap for the Gibson–Bruck next reaction method.
//! Gibson & Bruck (2000), DOI 10.1021/jp993732q.
#[doc(hidden)]
#[derive(Default)]
pub struct Clocks {
    heap: Vec<usize>,
    positions: Vec<usize>,
    times: Vec<f64>,
    rates: Vec<f64>,
}

impl Clocks {
    pub fn reset(&mut self, channels: usize) {
        self.heap.clear();
        self.heap.extend(0..channels);
        self.positions.clone_from(&self.heap);
        self.times.clear();
        self.times.resize(channels, f64::INFINITY);
        self.rates.clear();
        self.rates.resize(channels, 0.0);
    }

    pub fn next(&self) -> Option<(f64, usize)> {
        let i = *self.heap.first()?;
        self.times[i].is_finite().then_some((self.times[i], i))
    }

    pub fn needs_draw(&self, i: usize, rate: f64, fired: bool) -> bool {
        rate > 0.0 && (fired || self.rates[i] == 0.0)
    }

    /// Return false if the requested positive clock cannot advance in f64.
    /// Supply a fresh unit exponential for a fired or newly enabled channel.
    pub fn update(&mut self, i: usize, rate: f64, now: f64, fresh: Option<f64>) -> bool {
        if fresh.is_none() && rate == self.rates[i] {
            return true;
        }
        let old = self.times[i];
        let at = if rate == 0.0 {
            f64::INFINITY
        } else if let Some(unit) = fresh {
            now + unit / rate
        } else {
            now + (self.rates[i] / rate) * (old - now)
        };
        if rate > 0.0 && (!at.is_finite() || at <= now) {
            return false;
        }
        self.times[i] = at;
        self.rates[i] = rate;
        let mut pos = self.positions[i];
        if at < old {
            while pos > 0 {
                let parent = (pos - 1) / 2;
                if !self.less(pos, parent) {
                    break;
                }
                self.swap(pos, parent);
                pos = parent;
            }
        } else {
            loop {
                let left = pos * 2 + 1;
                if left >= self.heap.len() {
                    break;
                }
                let right = left + 1;
                let child = if right < self.heap.len() && self.less(right, left) {
                    right
                } else {
                    left
                };
                if !self.less(child, pos) {
                    break;
                }
                self.swap(pos, child);
                pos = child;
            }
        }
        true
    }

    fn less(&self, a: usize, b: usize) -> bool {
        let a = self.heap[a];
        let b = self.heap[b];
        (self.times[a], a) < (self.times[b], b)
    }

    fn swap(&mut self, a: usize, b: usize) {
        self.heap.swap(a, b);
        self.positions[self.heap[a]] = a;
        self.positions[self.heap[b]] = b;
    }
}
