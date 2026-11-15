# Allelescope user guide

Allelescope is a browser-based viewer for inspecting called variants and the
underlying alignment reads.  The page has two areas: the left pane is where you
choose a site and refine the reads; the right pane is the IGV alignment view.

## Open a browser instance

Start the server from the directory that contains one or more generated browser
instances:

```bash
start_server --data-dir /path/to/Browsers --port 8000
```

Open the instance in a browser.  For an instance in
`/path/to/Browsers/sample-01`, either of these addresses works:

```
http://localhost:8000/sample-01/
http://localhost:8000/sample-01/interactive_variants.html
```

If the server is on a remote machine, create an SSH tunnel first (see the
[README](README.md#accessing-from-a-local-machine-via-ssh-tunnel)).

## Find and inspect a called variant

1. In **Variants**, select a chromosome tab.  The number in parentheses is the
   number of called variants on that sequence.  For genomes with many
   sequences (segmented viruses, assemblies with many contigs), choose how the
   list is shown under **Sequence list** in the Variants gear menu:
   * **Automatic** (default): tabs across the top while they all fit in one
     row, otherwise tabs down the side;
   * **Tabs across the top**: one row that scrolls sideways;
   * **Tabs down the side**: a column beside the variant table, with a
     **Find sequence** box once there are 9 or more sequences; long names are
     shortened, and pointing at a tab shows the full name;
   * **Drop-down list**: one compact list.

   Your browser remembers the choice.
2. Click a row in the `Pos / Ref / Alt / Info` table.  The row is highlighted
   and the alignment view moves to that site.
3. Use **Flank (bp)** to choose how much sequence is shown on either side of a
   selected site (the default is 100 bp).  Click the row again after changing
   it to recenter with the new flank.
4. Inspect the reads in the right-hand IGV view.  Use its usual pan and zoom
   controls to look more closely at the local alignment.

The `Info` cell is editable, so it can be used for a temporary annotation while
reviewing a variant.  These edits are in the page only; they do not update the
VCF or any server-side file.

### Very deep sites: subsampling

The first time you open a site, Allelescope prepares its reads (the reads
crossing the site, their mates and split partners). Before doing this it
counts the read pairs crossing the site. If there are more than 2,000, it asks:

* **Subsample** shows the number of read pairs in the box (2,000 is suggested;
  you can type another number). Pairs are chosen at random, but the choice is
  fixed: a read always comes with its mate and split segments, and the same
  reads are shown each time the site is opened.
* **Use all reads** shows every read pair. At very high depth (for example a
  viral sample at 100,000x or more) this takes noticeably longer to prepare and
  to draw, and the view is crowded.

Tick **Do the same for other deep sites in this session** to stop being asked
again until the page is reloaded.

The allele counts in **Alleles at _chromosome:position_** always use every
read, whichever you choose. When a site is subsampled the panel says so, for
example *Subsampled: showing 2,000 of 195,351 read pairs*, with a **Show all
reads** button. A site shown in full has a **Subsample to 2,000** button.
Your browser remembers the choice for each site of each instance.

Filters, **Retain**/**Filter** of selected reads and **Export** work on the
reads being shown. The subsampled and full versions of a site are kept
separately on the server, so switching back to one you have already opened is
quick.

IGV has its own display limit: in a window with very many reads it draws only a
sample of them (by default 500 reads per 100 bases), though its coverage track
counts all of them. Subsampling here means the reads you select, filter and
export are the same set the server prepared.

The **Spanning Track** option shows or hides the `Spanning Reads Only` alignment
track.  **View as Pairs** changes that track to IGV's paired-read display.
The slider above the viewer changes the height of the spanning-read track.

## Move to a position not in the VCF

Use **Jump to position** at the top of the left pane.  Enter one or more
one-based positions, separated by spaces, commas, or semicolons, then select
**Add** (or press Enter).  The positions are added to the currently selected
chromosome tab as `Manual` rows and are immediately opened in the viewer.

Choose **Clear** to remove the manually added rows from the current chromosome.
Manual positions are useful for checking a suspected site that has not been
called.  They are not written back to the VCF.

### What is remembered

Your browser remembers the changes you make to the variant list: the positions
you add and the rows you remove with **X**.  It does not store the list itself.
Each instance applies those changes to its own VCF variants, so on one server:

* a position added by hand appears in every instance that has that sequence;
* a removed row is hidden in every instance that has a variant at that
  sequence and position;
* each instance otherwise keeps showing its own calls.

**Clear Persistent Variants** in the header menu forgets all of these changes.
Loading a saved variant list (Variants gear menu) replaces them with the
difference between that list and the instance's own variants.

## Read filtering

Changing a filter reloads the spanning-read track around the active variant.
The **Active Filters** card records the currently non-default settings.  The
gear menu on the Filters card can reset filters or save, load, and delete named
filter presets.  Presets live in the browser's local storage, so they are
available only in that browser profile.

The filter groups are:

| Group | Use it to |
| --- | --- |
| **Filter by Value** | Limit reads by soft-clip length, edit distance, or mapping quality. |
| **BAQ Filters** | Limit reads by their base alignment quality at the selected site. |
| **Include (Require)** | Keep reads that have selected SAM flags or `XO` tag bits. |
| **Exclude (Forbid)** | Remove reads that have selected SAM flags or `XO` tag bits. |
| **Strand** | Keep reads from both strands, the plus (forward) strand, or the minus (reverse) strand. |

For the flag and `XO` tag checkboxes, **Include** and **Exclude** are mutually
exclusive for the same property: selecting one clears the opposite setting.
Use Include when that property must be present; use Exclude when it must be
absent.  The `XO` choices (`overlap`, `mates`, `split mates`, `split reads`,
`proper`, and `improper`) apply to Allelescope's `XO` read tag, when present.
`proper` means the aligner flagged the pair as proper and the mates face each
other, whichever of read 1 and read 2 is on the forward strand (`F1R2` or
`R1F2`); every other paired read is `improper`.

### Numeric limits

The numeric boxes are empty by default, and an empty box means no limit, so
nothing is filtered until you type a value.  A value is applied as typed,
including 0: a maximum edit distance of 0 keeps only reads that match the
reference exactly.

| Box | Keeps reads with |
| --- | --- |
| **Soft Clip Length** Min / Max | at least / at most that many soft-clipped bases (both ends added together) |
| **Edit Distance** Min / Max | an `NM` edit distance of at least / at most that value |
| **Mapping Quality** Min | a mapping quality of at least that value |
| **BAQ value** Min / Max | a base alignment quality at the selected site of at least / at most that value |

BAQ is read from the `BQ` tag that `samtools calmd -r` adds while a site is
prepared.  It is evaluated on the reads that cover the site, and a read pair is
kept or removed as a whole, so mates follow the read that covers the site.  For
a read with a deletion at the site, the median of the neighbouring bases is
used.

The **Soft clips** tick boxes are separate from the Soft Clip Length boxes:
under Include the read must have at least one soft-clipped base, under Exclude
it must have none.

### Active Filters

Each badge states the condition that is applied.  `require ...` means only
reads with that property are shown (Include ticks and minimums, for example
`require MAPQ ≥ 20`); `exclude ...` means reads with that property are hidden
(Exclude ticks and maximums, for example `exclude duplicates` or
`exclude soft clip > 5`).  Click the × on a badge to remove that filter.

## Filter by observed allele at a site

When you select a variant row or add a manual position, Allelescope extracts
the observed alleles and shows an **Alleles at _chromosome:position_** panel
below the variant list.  Each row reports an allele type (for example `REF`,
`SNP`, `DEL`, or `INS`), its supporting read-pair count, and its percentage.
For a variant from the VCF, the row that matches the record's ALT allele is
marked `ALT`.

A VCF records an insertion or deletion at the base before the change, where
every read still shows the reference.  For those records Allelescope compares
the reads one base further on (the first deleted base, or the base after the
insertion point) and says so in the panel; the view stays centred on the VCF
position.

1. Choose which allele rows are visible: **All**, the **Top** _N_, or those at
   least a chosen percentage.
2. Tick the alleles of interest, or use **All**/**None** to set the checkboxes.
3. Select **Filter** to retain reads associated with the selected allele(s).

This allele filter is combined with the other active filters and is listed
under **Active Filters** (for example `require allele DEL 7`).  It belongs to
the site it was set at: selecting another site clears it.  To clear it
yourself, remove its badge, use **Reset filters**, or choose **Filter** with no
allele ticked.

## Select, compare, and export reads

To select a read from the alignment view, hold **Shift** and click the read.
Its name appears in **Read Selection** and the selected reads are highlighted.
You can also paste a read name into the field and choose **Add**.

For selected reads:

| Action | Result |
| --- | --- |
| Click a read name | Searches for that read in IGV. |
| **Retain** | Restricts the spanning-read track to the selected reads. |
| **Filter** | Removes the selected reads from the spanning-read track. |
| **Export** | Downloads the selected reads after the active filters are applied. |
| **Clear List** | Clears the selection and restores the unselected read set. |

**Retain** and **Filter** apply to the site you are viewing.  When you select
another site the full read set is shown again; the selection is kept, so you
can apply it there too.

> **Current limitation:** the **Clear Highlights** button is displayed, but its
> action is not implemented in the current viewer.  Use **Clear List** when you
> need to clear the visual highlighting as well as the selection.

## Make more room for the alignment view

Drag the divider between the panes to resize the left pane.  Select the arrow
button at its lower edge to collapse or reopen it; the browser remembers that
collapsed state locally.  The viewer remains usable with the left pane hidden
once a site has been opened.

## Practical review sequence

For a typical manual review, select the called variant, inspect both strands,
then apply one filter at a time (for example, remove duplicates or require
proper pairs).  Select suspicious or supporting reads with Shift-click,
compare them using **Retain** or **Filter**, and export the final selection if
needed.  Reset filters before moving on when you want to return to the full
read set.

## TODO
reads are colored by strand