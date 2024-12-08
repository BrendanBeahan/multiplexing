import argparse
import pandas as pd
import numpy as np
import pysam
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import matplotlib.pyplot as plt

# Step 1: Map reference names to barcodes
ref2bc = {
    'cc6m_2244_T7_ecorv': 1,
    'cc6m_2459_T7_ecorv': 2,
    'cc6m_2595_T7_ecorv': 3,
    'cc6m_2709_T7_ecorv': 4
}

def main(bam_file, seqtagger_output, output_file=None):

    # Step 2: Extract true barcodes from BAM file
    read2bc = {}
    sam = pysam.AlignmentFile(bam_file)

    for a in sam:
        if a.is_secondary or a.mapq < 20 or a.is_supplementary:
            continue
        r, ref = a.qname, a.reference_name
        if ref not in ref2bc:
            continue
        read2bc[r] = ref2bc[ref]

    print(f"Total reads with true barcodes: {len(read2bc)}")

    # Step 3: Load SeqTagger output and prepare data
    df = pd.read_csv(seqtagger_output, sep='\t')

    # Extract predicted barcodes
    predicted_barcodes = df['barcode'].astype(int)
    read_ids = df['read_id']

    # Map read IDs to true barcodes
    true_barcodes = [read2bc.get(rid, 0) for rid in read_ids]

    # Filter out reads without true barcodes
    valid_idx = [i for i, bc in enumerate(true_barcodes) if bc > 0]
    y_true = [true_barcodes[i] for i in valid_idx]
    y_pred = [predicted_barcodes[i] for i in valid_idx]

    # Step 4: Generate Confusion Matrix
    labels = sorted(ref2bc.values())
    cm = confusion_matrix(y_true, y_pred, labels=labels)

    # Step 5: Display Confusion Matrix
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=labels)
    disp.plot(cmap='Blues', values_format='d')
    plt.title('Confusion Matrix: SeqTagger')
    plt.xlabel('Predicted Barcode')
    plt.ylabel('True Barcode')

    # Save Confusion Matrix to a file if output_file is provided
    if output_file:
        plt.savefig(output_file)
        print(f"Confusion matrix saved to {output_file}")

    # Display the plot
    plt.show()

    # Optional: Print Accuracy
    accuracy = np.trace(cm) / np.sum(cm)
    print(f"Overall Accuracy: {accuracy:.2%}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate confusion matrix for SeqTagger output.')
    parser.add_argument('-b', '--bam', required=True, help='Path to the BAM file.')
    parser.add_argument('-s', '--seqtagger', required=True, help='Path to the SeqTagger output TSV file.')
    parser.add_argument('-o', '--output', required=False, help='Path to save the confusion matrix image (optional).')
    args = parser.parse_args()

    main(args.bam, args.seqtagger, args.output)
