#include "cigar.h"

uint32_t cigarToInt(char operation, uint32_t len) {
	uint32_t opCode;
	// Convert CIGAR operation character to an operation code
	switch (operation) {
	case 'M': opCode = 0x0; break; // Match
	case 'I': opCode = 0x1; break; // Insertion
	case 'D': opCode = 0x2; break; // Deletion
	case '=': opCode = 0x7; break; // Sequence match
	case 'X': opCode = 0x8; break; // Mismatch
		// Add cases for other SAM specification operation codes as needed
	default: opCode = 0xF; break; // Unknown operation
	}
	// Combine operation length and code into a single uint32_t value
	return (len << 4) | opCode; // Shift length left by 4 bits, then combine with opCode
}

void intToCigar(CigarUnit cigar, char& operation, uint32_t& len)
{
	uint32_t opCode = cigar & 0xF; // Extract the lower 4 bits as the operation code
	len = cigar >> 4; // Extract the length by shifting right by 4 bits

	// Convert operation code back to a CIGAR operation character
	switch (opCode) {
	case 0x0: operation = 'M'; break; // Match
	case 0x1: operation = 'I'; break; // Insertion
	case 0x2: operation = 'D'; break; // Deletion
	case 0x7: operation = '='; break; // Sequence match
	case 0x8: operation = 'X'; break; // Mismatch
		// Add cases for other SAM specification operation codes as needed
	default: operation = '?'; break; // Unknown operation
	}
}

Score_t caculateMatchScore(const char* match, uint_t length)
{
	Score_t score = 0;
	
	for (size_t i = 0; i < length; ++i) {
		int id = ScoreChar2Idx[match[i]];
		if (id == 0 || id == 3) {
			score += MATCH_SCORE_AT;
		}
		else if (id == 1 || id == 2) {
			score += MATCH_SCORE_CG;
		}
		else {
			score += MATCH_SCORE_N;
		}
	}
	return score;
}
