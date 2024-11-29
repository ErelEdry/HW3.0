import math

# Example Sequence Corpus
seq_corpus = {
    1: 'GCTTA_tGC:GXATC  CGTAGACffx:TYAGgytACGTMA',
    2: 'AGG. Ddfe::wscv',
    3: 'cl_yuCATGATGCGTACCAGGCTqwAGCATGCGTbbAGCTAxzvGCATGAC'
}


# Helper function
def deep_copy_dict(original_dict):
    """
    This function copies the dictionary.

    :return: A new dictionary that is a deep copy of the original.
    """
    copied_dict = {}  # a empty dictionary to hold the copied values.

    for key, value in original_dict.items():
        if isinstance(value, dict):  # If the value is a dictionary
            copied_dict[key] = deep_copy_dict(value)  # Recursively copy nested dictionaries
        elif isinstance(value, list):  # If the value is a list
            copied_dict[key] = value.copy()  # Perform a shallow copy of the list
        else:
            copied_dict[key] = value  # For non-dictionary, non-list values, copy directly

    return copied_dict


def codon_translator(codon):
    RNA_to_protien = {'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
                      'UAU': 'Y', 'UAC': 'Y',
                      'UGU': 'C', 'UGC': 'C', 'UGG': 'W', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L', 'CCU': 'P',
                      'CCC': 'P', 'CCA': 'P',
                      'CCG': 'P', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R',
                      'CGG': 'R', 'AUU': 'I',
                      'AUC': 'I', 'AUA': 'I', 'AUG': 'M', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAU': 'N',
                      'AAC': 'N', 'AAA': 'K',
                      'AAG': 'K', 'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V',
                      'GUG': 'V', 'GCU': 'A',
                      'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGU': 'G',
                      'GGC': 'G', 'GGA': 'G',
                      'GGG': 'G'}
    decoding = RNA_to_protien.get(codon, "stop")
    return decoding


####### Part A #########
def remove_punctuation(seq):
    """
    this function remove punctuation from the sequence
    :return: str
    """
    punctuation = ":-_."
    str = ""
    i = 0
    while i < len(seq):
        if seq[i] not in punctuation:  # if the character not in the punctuation
            str += seq[i]  # adding the result to the new string
        else:
            str += " "  # if is a punctuation replace with space
        i += 1
    return str


def remove_spaces(seq):
    """
    this function remove spaces
    :return: str
    """
    str = ""
    i = 0
    while i < len(seq):
        if seq[i] != " ":  # if the character not is not sapce
            str += seq[i]  # adding the result to the new string
        i += 1
    return str


def remove_invalid_letters(seq):
    """
    this function remove invalid letters
    :return: str
    """
    valid_letters = "ATGCatgc"
    str = ""
    i = 0
    while i < len(seq):
        if seq[i] not in valid_letters:
            str += " "  # adding the result to the new string
        else:
            str += seq[i]
        i += 1

    return str


def capitalize_letters(seq):
    """
    this function making sure that the letters are capitalized
    :return: capitalize_str
    """
    capitalize_str = ""  # temp str
    for i in range(len(seq)):
        if seq[i].islower():  # checking if the character is lowercase
            capitalize_str += seq[i].upper()  # making it upper and inserting it to the temp str
        else:
            capitalize_str += seq[i]  # adding the good character to the temp str

    return capitalize_str


def proofreading(seq):
    """
    this function getting a string just with ATCG sequence and returning list of strings as third of functional DNA sequence
    :return:func_seq
    """
    valid_characters = ["A", "T", "C", "G"]
    stop_condos = ["TAA", "TAG", "TGA"]
    func_seq = []
    atg_index = -1
    for i in range(len(seq)):  # Checks if there are forbidden letters in the sequence
        if seq[i] not in valid_characters:
            return []
    for i in range(len(seq) - 2):  # we need to start with ATG
        if seq[i:i + 3] == "ATG":
            atg_index = i
            break
    seq = seq[atg_index:]
    if seq[:3] == "ATG":  # validate that the str starts with ATG
        if len(seq) % 3 != 0:  # if the seq not dividing by 3
            seq = seq[:-(len(seq) % 3)] + "TAA"
            for i in range(0, len(seq), 3):
                func_seq.append(seq[i:i + 3])  # Slices a relevant three and adds as an element in the list
        else:
            if seq[-3:] not in stop_condos:  # checking if the last 3 letters is from the stop condos
                seq = seq[:-3] + "TAA"
                for i in range(0, len(seq), 3):
                    func_seq.append(seq[i:i + 3])  # Slices a relevant three and adds as an element in the list
            else:
                for i in range(0, len(seq), 3):
                    func_seq.append(seq[i:i + 3])  # Slices a relevant three and adds as an element in the list

    return func_seq


def transcript(func_seq):
    """
    this function translate dna to rna
    :return: func_seq
    """
    for i in range(len(func_seq)):
        new_str = ""  # temp string
        for j in range(len(func_seq[i])):
            if func_seq[i][j] == "T":  # if we see T
                new_str += "U"  # insert a U to the temp string instead
            else:
                new_str += func_seq[i][j]
        func_seq[i] = new_str  # adding the transcripted str to the list

    return func_seq


def translate(rna_seq):
    """
    this function taking sequence to amino acid
    :return: rna_list
    """
    rna_list = []
    for i in range(len(rna_seq)):
        if codon_translator(rna_seq[i]) != 'stop':
            rna_list.append(codon_translator(rna_seq[i]))#adding to the list
        if codon_translator(rna_seq[i]) == 'stop': #if we getting a stop, we stopping the function
            break

    return rna_list


def preprocessing(seq):
    """using all the functions by the order that we asked for"""
    seq = remove_punctuation(seq)
    seq = remove_invalid_letters(seq)
    seq = remove_spaces(seq)
    seq = capitalize_letters(seq)
    return translate(transcript(func_seq=proofreading(seq)))

###### Part B #######
def get_sequences_data(sequence_corpus):
    """
    this function gets sequences data from the sequence corpus
    :return: sequences_data
    """
    sequences_data = {}
    pro_corpus = {key: preprocessing(sequence_corpus[key]) for key in
                  sequence_corpus}  # copy to keep the dic like at the start
    for key in pro_corpus:
        sequences_data[key] = len(pro_corpus[key])  # the value is the length of the amino acid sequences

    return sequences_data


def create_inverted_index(sequence_corpus):
    """
    this function creates inverted index
    :return: inverted_index
    """
    inverted_index = {}
    for key in sequence_corpus:
        sequence_corpus[key] = preprocessing(sequence_corpus[key])  # proccecing the dna
    for key, value_list in sequence_corpus.items():  # Going through both the value and the key at the same time
        for char in value_list:
            if char not in inverted_index:  # Add an empty dictionary when the rest of the letter does not exist
                inverted_index[char] = {}
            if key in inverted_index[char]:  # If the key exists we will add the number of its occurrences accordingly
                inverted_index[char][key] += 1
            else:
                inverted_index[char][key] = 1

    return inverted_index


def add_to_data(inverted_index, sequences_data, seq_id, seq):
    """
    this function adding sequence to the data
    :return: sequences_data,inverted_index
    """
    pro_seq = preprocessing(seq)  # Process the sequence
    sequences_data[seq_id] = len(pro_seq)  # Adding the seq_id as a key, and the len of processed sequence
    for i in pro_seq:  # Iterate over each character in the processed sequence
        if i not in inverted_index:  # Check if the character is not exist
            inverted_index[i] = {}  # Add an empty dictionary for this character

        if seq_id in inverted_index[i]:  # Check if seq_id exists for this character
            inverted_index[i][seq_id] += 1  # adding 1
        else:
            inverted_index[i][seq_id] = 1  # If seq_id not exist,put 1

    return sequences_data, inverted_index


def remove_from_data(inverted_index, sequences_data, seq_id):
    """
    this function removes sequence from the data
    :return: inverted_index,sequences_data
    """
    if seq_id in sequences_data:  # checikng if the seq_id exists
        sequences_data.pop(seq_id)  # using pop function to remove it
    for key, value in inverted_index.items():
        if seq_id in value:  # checikng if the seq_id exists
            value.pop(seq_id)  # using pop function to remove it

    return inverted_index, sequences_data


def preprocess_query(query):
    """
    this function taking query and preprocesses it to a tuple
    :return: query_tuple
    """
    query_tuple = tuple(query.replace(",", ""))  # creating the tuple and deleting all the commas

    return query_tuple


###### Part C #######
def calculate_aaf_isf(amino_acid, seq_id, inverted_index, sequences_data):
    """
    The function calculates the AAF and ISF values and calculates ISF*AAF
    """
    shallow_inverted_index = deep_copy_dict(inverted_index)  # Perform a deep copy of the inverted_index dictionary to prevent modifying the original dictionary
    shallow_sequences_data = deep_copy_dict(sequences_data)  # Perform a deep copy of the sequences_data dictionary to prevent modifying the original dictionary
    if shallow_sequences_data.get(seq_id, 0) != 0:
        AAF = (math.log(len(shallow_sequences_data) / len(shallow_inverted_index.get(amino_acid, {})),2))  # calculation of AAF
        ISF = shallow_inverted_index.get(amino_acid, {}).get(seq_id, 0) / shallow_sequences_data.get(seq_id,1)  # calculation of ISF
        ans = AAF * ISF
        return round(ans, 3)

def get_scores_of_relevance_sequences(query, inverted_index, sequences_data):
    """
        this function getting a tuple and returning a dictionary with scores of aaf-isf between the amino acids to the sequences
    :return: scores_seq
    """
    scores_seq = {}  # Create an empty dictionary to store the final results

    shallow_inverted_index = deep_copy_dict(inverted_index)  # Perform a deep copy of the inverted_index dictionary to prevent modifying the original dictionary
    shallow_sequences_data = deep_copy_dict(sequences_data)  # Perform a deep copy of the sequences_data dictionary to prevent modifying the original dictionary

    for amino_acid in query:
        if amino_acid in shallow_inverted_index.keys():  # If the amino acid exists in the shallow_inverted_index dictionary
            for seq_id in shallow_inverted_index[amino_acid]:
                if seq_id not in scores_seq:  # If the sequence ID is not already in the scores_seq dictionary, put 0
                    scores_seq[seq_id] = 0
                scores_seq[seq_id] = round(
                    scores_seq[seq_id] + calculate_aaf_isf(amino_acid, seq_id, inverted_index, shallow_sequences_data),
                    3)  # Calculate the AAF-ISF score and add it to the scores_seq dictionary
    return scores_seq

sequences_data = get_sequences_data(seq_corpus)
inverted_index = create_inverted_index(seq_corpus)
print(get_scores_of_relevance_sequences(('R', 'M', 'S'),
inverted_index, sequences_data))
###### Part D #######
def menu(sequence_corpus):
    """
    Displays the menu and allows the user to choose an option.
    Handles different actions based on the user's choice.
    """
    choice = input('Choose an option from the menu:\n\t(1) Insert a query.\n\t(2) Add sequence to sequence_corpus.\n\t(3) Calculate AAF-ISF Score for an amino acid in a sequence.\n\t(4) Delete a sequence from the sequence_corpus.\n\t(5) Exit.\nYour choice: ')


    if choice == '1':  # Insert a query
        query_choice = input('Choose the type of results you would like to retrieve:\n\t(A) All relevant sequences.\n\t(B) The most relevant sequence.\n\t(C) Back to the main menu.\nYour choice: ')
        query = input('Write your query here:')
        if query_choice == 'A':  # Show all relevant sequences
            print(get_scores_of_relevance_sequences(preprocess_query(query), create_inverted_index(sequence_corpus),get_sequences_data(sequence_corpus)))
        elif query_choice == 'B':  # Show the most relevant sequence
            highest_score_seq_id=max(calculate_aaf_isf(), key=calculate_aaf_isf().get)
            highest_score=max(calculate_aaf_isf().values())
            print(f'The most relevant sequence is {highest_score_seq_id} with a score of {highest_score}')
        elif query_choice == 'C':  # Back to main menu
            menu(sequence_corpus)
        else:
            print('Invalid choice. Please select a valid option.')
            menu(sequence_corpus)
    elif choice == '2':  # Add sequence to corpus
        while True:
            seq_id = int(input("Insert the sequence ID:"))
            if seq_id not in sequence_corpus.keys():  # Check if sequence is already present
                new_seq=input(f'Insert the sequence itself:')
                add_to_data( create_inverted_index(sequence_corpus),get_sequences_data(sequence_corpus), seq_id, new_seq)
                print(f'Sequence {seq_id} was successfully added!')
                break
            else:
                print(f'The sequence ID {seq_id} is already in sequence_corpus.')
    elif choice == '3':  # Calculate AAF-ISF Score for an amino acid
        while True:
            seq_id = int(input("Insert the sequence ID:"))
            key_list=list(get_sequences_data(seq_corpus).keys())
            if int(seq_id) not in key_list :  # Check if sequence exists
                print(f'The sequence ID {seq_id} is not in sequence_corpus.')
            else:
                amino_acid_after_preprocces = input("Insert the amino acid:")
                while True:
                    if amino_acid_after_preprocces not in preprocessing(sequence_corpus[seq_id]):  # Check if amino acid is valid
                        print(f'The amino acid {amino_acid_after_preprocces} is not in sequence_corpus.')
                    else:
                        aaf_isf = calculate_aaf_isf(amino_acid_after_preprocces, seq_id,  create_inverted_index(sequence_corpus),get_sequences_data(sequence_corpus))
                        print(
                            f'AAF-ISF of the amino acid {amino_acid_after_preprocces} in sequence {seq_id} is: {aaf_isf}')
                        break
    elif choice == '4':  # Delete a sequence from corpus
        seq_id = int(input("Insert the sequence ID:"))
        while True:
            if seq_id not in get_sequences_data(seq_corpus).keys():  # Check if sequence exists
                print(f'The sequence ID {seq_id} is not in sequence_corpus.')
            else:
                remove_from_data( create_inverted_index(sequence_corpus),get_sequences_data(sequence_corpus), seq_id)  # Remove sequence
                print(f'Sequence {seq_id} was successfully deleted.')
                break
    elif choice == '5':  # Exit
        return  # Finish the function
    elif choice not in ["1", "2", "3", "4", "5"]: #Option for invalid strings
        print( 'Invalid choice. Please select a valid option.')
        menu(sequence_corpus)




