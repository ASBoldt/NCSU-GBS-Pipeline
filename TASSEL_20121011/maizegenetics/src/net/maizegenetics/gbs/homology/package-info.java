/**
 * Homology package has methods for identifying homology between sequence tags.
 * <p>
 * ParseBarcodeRead, Barcode, and ReadBarcaodeResult support the processing of barcodes.
 * <p>
 * TagMatchFinder & ReadBLASTer are both BLAT like matching algorithms (Ed - I am
 * not sure how they are different).  Smith-Waterman is the basic alignment algorithm,
 * however, some of the biojava implementations may be better.
 * <p>
 * TagToAlleles is a data structure for storing these homology relationships.
 *
 *
 */

package net.maizegenetics.gbs.homology;
