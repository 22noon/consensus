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

Open the instance's `interactive_variants.html` page in a browser.  For an
instance in `/path/to/Browsers/sample-01`, the address is:

```
http://localhost:8000/sample-01/interactive_variants.html
```

If the server is on a remote machine, create an SSH tunnel first (see the
[README](README.md#accessing-from-a-local-machine-via-ssh-tunnel)).

## Find and inspect a called variant

1. In **Variants**, select a chromosome tab.  The number in parentheses is the
   number of called variants on that sequence.
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

## Read filtering

Changing a filter reloads the spanning-read track around the active variant.
The **Active Filters** card records the currently non-default settings.  The
gear menu on the Filters card can reset filters or save, load, and delete named
filter presets.  Presets live in the browser's local storage, so they are
available only in that browser profile.

The filter groups are:

| Group | Use it to |
| --- | --- |
| **Filter by Value** | Limit reads by soft-clip length or edit distance. |
| **Include (Require)** | Keep reads that have selected SAM flags or `XO` tag bits. |
| **Exclude (Forbid)** | Remove reads that have selected SAM flags or `XO` tag bits. |
| **Strand** | Keep reads from both strands, the plus strand, or the minus strand. |

For the flag and `XO` tag checkboxes, **Include** and **Exclude** are mutually
exclusive for the same property: selecting one clears the opposite setting.
Use Include when that property must be present; use Exclude when it must be
absent.  The `XO` choices (`overlap`, `mates`, `split mates`, `split reads`,
`proper`, and `improper`) apply to Allelescope's `XO` read tag, when present.

Soft-clip filtering uses the inclusive minimum field to retain reads with at
least that many clipped bases and the maximum field to retain reads with no
more than that many.  Edit Distance currently acts as an upper limit despite
the UI's `Min` label: a non-zero value retains reads whose `NM` edit-distance
tag is at most that value.

> **Current limitation:** the visible **Mapping Quality** and **BAQ value**
> fields are not applied by the current server filtering code.  Do not rely on
> them to change the displayed reads.

## Filter by observed allele at a site

When you add a manual position, Allelescope extracts the observed alleles and
shows an **Alleles at _chromosome:position_** panel.  Each row reports an allele
type (for example `REF`, `SNP`, `DEL`, or `INS`), its supporting read count, and
its percentage of the extracted reads.

1. Choose which allele rows are visible: **All**, the **Top** _N_, or those at
   least a chosen percentage.
2. Tick the alleles of interest, or use **All**/**None** to set the checkboxes.
3. Select **Filter** to retain reads associated with the selected allele(s).

This allele filter is combined with the other active filters.  If no allele is
selected, the viewer leaves the current read set unchanged and reports that
nothing was selected.

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
