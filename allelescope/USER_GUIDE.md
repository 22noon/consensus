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

The variant table has these columns:

| Column | Shows |
| --- | --- |
| status icon | Your review status: empty circle (not reviewed), Confirmed, Doubtful or Artefact.  Click it to change. |
| **Pos**, **Ref**, **Alt** | The VCF record.  Long alleles are shortened; hover to see them in full. |
| **AF** | The allele fraction from the VCF (`INFO/AF`, or the sample's `AF`, `VAF` or `AD`).  If the VCF has none, the fraction seen in the reads appears *in italics* once you have opened the site. |
| **Note** | Your note, or the VCF `INFO` summary until you write one.  Double-click to edit; Enter keeps it, Esc cancels. |
| × | Removes the row (it appears when you point at the row).  A message offers **Undo**. |

Click **Pos**, **AF** or the status icon in the header to sort; click again to
reverse.  Above the table, **Search** finds positions, alleles and notes, and
the list next to it shows only variants with one status (for example
**Not reviewed**).

### Reviewing a list of variants

Statuses and notes are kept in your browser, separately for each instance,
and survive reloading the page.  The count above the table shows how many
variants have a status.  The keyboard makes a long list quick to go through
(press **?** or the keyboard button for this list; shortcuts are off while you
type in a box):

| Key | Does |
| --- | --- |
| ↓ or J, ↑ or K | Next / previous variant |
| N | Next variant without a status, across all sequences |
| C, D, X | Mark Confirmed, Doubtful, Artefact (press again to clear) |
| U | Clear the status |
| E | Edit the note |
| 1 – 9 | Show only reads with that allele of the allele panel |
| 0 | Show reads with any allele |
| / | Search the variant list |

The download button above the table saves the review:

* **Review table (.tsv)**: every variant with its status, note, AF and, for
  the sites opened in this session, the alleles seen with their strand counts.
* **VCF with review (.vcf)**: the instance's VCF with `INFO/REVIEW`
  (`confirmed`, `doubtful` or `artefact`) and `INFO/REVIEW_NOTE` added to the
  reviewed records.  Notes are percent-encoded (`%20` for a space, `%3B` for
  `;`).  Positions you added by hand are not in the VCF; they are in the .tsv.

### Sharing a view

The address in the browser records what you are looking at: the variant, the
filters, the allele filter and the subsample size, for example
`.../demo/#v=demo_seg1:2450&on=retain-duplicates&xa=MC`.  Send the address (or
use the link button in the allele panel to copy it) and whoever opens it on
the same server sees the same view.  A position that is not in their list is
added to it.  The browser's **Back** and **Forward** buttons move between the
variants you have opened.

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

Filters, **Show only these**/**Hide these** for selected reads and **Export** work on the
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
The filters in use are listed at the top of the Filters card, and its header
shows how many there are.  The gear menu on the Filters card can reset
filters, download the settings as a file, or save, load, and delete named
filter presets.  Presets live in the browser's local storage, so they are
available only in that browser profile.

The filter groups are:

| Group | Use it to |
| --- | --- |
| **Filter by Value** | Limit reads by soft-clip length, edit distance, or mapping quality. |
| **BAQ Filters** | Limit reads by their base alignment quality at the selected site. |
| **Flags and read types** | Require or exclude SAM flags, Allelescope's `XO` read types and pair orientations. |
| **Strand** | Keep reads from both strands, the plus (forward) strand, or the minus (reverse) strand. |

Under **Flags and read types** each property has three settings: **any** (no
filter), **require** (only reads with it are shown) and **exclude** (reads with
it are hidden).  Point at a name to see what it means.  The `XO` choices (`overlap`, `mates`, `split mates`, `split reads`,
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

### The list of filters in use

Each badge states the condition that is applied.  `require ...` means only
reads with that property are shown (required properties and minimums, for example
`require MAPQ ≥ 20`); `exclude ...` means reads with that property are hidden
(excluded properties and maximums, for example `exclude duplicates` or
`exclude soft clip > 5`).  Click the × on a badge to remove that filter.

## Filter by observed allele at a site

When you select a variant row or add a manual position, Allelescope extracts
the observed alleles and shows an **Alleles at _chromosome:position_** panel
below the variant list.  Each row reports an allele type (for example `REF`,
`SNP`, `DEL`, or `INS`), its supporting read-pair count, and its percentage.
For a variant from the VCF, the row that matches the record's ALT allele is
marked `ALT`.

The line under each allele gives the reads on the forward and reverse strand
(`15+ / 11−`), the mean base quality at the site (`BQ`, not shown for
deletions) and the mean mapping quality (`MQ`).  An allele marked
**⚠ strand** is mostly on one strand while the other reads at the site are
not (Fisher's exact test p < 0.01, with at most 15 % on one strand), which is
a common sign of an artefact.  The counts use one read per read pair: the
read that crosses the site.

### Read evidence

Under the allele list, **Read evidence** compares the reads behind the main
alleles: how often the site is within 10 bp of a read end, other edits per
read (`NM`, not counting the allele), multi-mapped reads (MAPQ 0, secondary,
bwa `XA`/`XS`), soft clips near the site, split reads, median insert size,
F1R2 pair orientation and duplicates.  Each non-reference allele is compared
with the reference allele (or the commonest other allele), and a warning such
as **⚠ read end**, **⚠ mismatches**, **⚠ multi-mapped**, **⚠ orientation** or
**⚠ duplicates** is shown when the difference is both large and significant
(at least 5 reads, Fisher's exact test p < 0.01 where counts are compared).
Point at a warning to see what it means; the thresholds are in
`static/js/evidence.js` (`EVIDENCE_RULES`) and in the user manual, section 9.3.

The variant table has a flag column: a tick for an opened site with no
warnings, or ⚠ and the number of warnings, including **AF ≠ VCF** when the
VCF's AF and the reads differ by 15 points or more.  The warnings are added to
the review table (`read_flags`) and the reviewed VCF (`INFO/READ_FLAGS`).

Click an allele to show only its reads, and click it again to show all reads
(or press its number, 1 to 9, and 0 to clear).  While a site is being
prepared for the first time the panel shows *Preparing reads at ...*.

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
Its name appears in the **Selected reads** bar above the alignment view and
the selected reads are highlighted.  You can also paste a read name into the
field and choose **Add** (or press Enter).

For selected reads:

| Action | Result |
| --- | --- |
| Click a read name | Searches for that read in IGV. |
| **Show only these** | Restricts the spanning-read track to the selected reads. |
| **Hide these** | Removes the selected reads from the spanning-read track. |
| **Export** | Downloads the selected reads after the active filters are applied. |
| **Clear List** | Clears the selection and restores the unselected read set (with **Undo**). |

**Show only these** and **Hide these** apply to the site you are viewing.  When you select
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

For a typical manual review, press **N** to open the first variant without a
status.  Check the allele panel: the strand split and qualities of the ALT
allele, and any **⚠ strand** mark.  Press the ALT allele's number to see only
its reads, then apply one filter at a time (for example, exclude duplicates or
require proper pairs).  Select suspicious or supporting reads with
Shift-click and compare them with **Show only these** or **Hide these**.  Mark
the variant with **C**, **D** or **X**, add a note with **E** if needed, and
press **N** for the next one.  Download the review when you are done.

## Stored data and settings

The gear menu at the top right of the page says what each action affects:

* **Server: cached filtered reads**: delete the filtered reads cached for this
  instance or for all instances.  They are rebuilt when needed; prepared sites
  are kept.
* **This browser**: forget the added and removed variants (all instances),
  forget the review statuses and notes of this instance, or reset everything
  Allelescope saved in this browser.  The first two can be undone straight
  away; the reset cannot.

A key above the alignment view shows the read colours: pink reads are on the
forward strand and blue reads on the reverse strand; mismatched bases are
green (A), blue (C), orange (G) and red (T), a purple bar is an insertion and
a black line a deletion.  IGV draws at most 500 reads per 100 bases; its
coverage bar counts all of them.
