//! `BinSource`: trait abstracting over the two readers a `BinScanner`
//! can run against — the in-memory [`crate::BinTable`] (built up at
//! `crispr-ots build` time, also used by integration tests) and the
//! mmap'd zero-copy [`crate::MmapDb`] (used by `crispr-ots discover`).
//!
//! All the per-scan work happens off the trait: the scanner pulls the
//! enzyme + bin width + a snapshot of bin offsets at the very start of
//! a scan, then never touches the trait again inside the hot loop. So
//! dynamic dispatch via `&dyn BinSource` is fine — the cost is
//! one virtual call per top-level `scan()` invocation.

use crispr_encoding::{Enzyme, Site};

use crate::packed::PackedEntry;

/// Common read interface for "something a scanner can scan."
pub trait BinSource: Sync {
    fn enzyme(&self) -> &Enzyme;
    fn bin_width(&self) -> u8;
    fn bin_key(&self, site: Site) -> u32;
    fn contig_name(&self, id: u32) -> Option<&str>;
    /// Iterate `(bin_key, packed_entries)` over **non-empty** bins. The
    /// scanner collects this once at the top of each `scan()` call.
    fn iter_bins<'a>(&'a self) -> Box<dyn Iterator<Item = (u32, &'a [PackedEntry])> + 'a>;
}

impl BinSource for crate::BinTable {
    fn enzyme(&self) -> &Enzyme {
        self.enzyme()
    }
    fn bin_width(&self) -> u8 {
        self.bin_width()
    }
    fn bin_key(&self, site: Site) -> u32 {
        self.bin_key(site)
    }
    fn contig_name(&self, id: u32) -> Option<&str> {
        self.contig_name(id)
    }
    fn iter_bins<'a>(&'a self) -> Box<dyn Iterator<Item = (u32, &'a [PackedEntry])> + 'a> {
        Box::new(self.iter())
    }
}

impl BinSource for crate::MmapDb {
    fn enzyme(&self) -> &Enzyme {
        self.enzyme()
    }
    fn bin_width(&self) -> u8 {
        self.bin_width()
    }
    fn bin_key(&self, site: Site) -> u32 {
        self.bin_key(site)
    }
    fn contig_name(&self, id: u32) -> Option<&str> {
        self.contig_name(id)
    }
    fn iter_bins<'a>(&'a self) -> Box<dyn Iterator<Item = (u32, &'a [PackedEntry])> + 'a> {
        Box::new(self.iter())
    }
}
