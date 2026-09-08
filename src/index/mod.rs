//! On-disk target index: the mmapped archive, its suffix array, and the views
//! the search reads it through.

mod archive;
pub(crate) mod sa;
pub(crate) mod store;
pub(crate) mod view;

pub use view::TargetView;
