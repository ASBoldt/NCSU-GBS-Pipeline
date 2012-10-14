package net.maizegenetics.gbs.pipeline;

/** Stores a single record from a TagCount file.
 * @author jvh39 */
public class TagCountRecord {
    long[] sequence;
    byte length;
    int count;

    TagCountRecord(byte tagLengthInLong){
        sequence = new long[tagLengthInLong];
        length = 0;
        count = 0;
    }
}
