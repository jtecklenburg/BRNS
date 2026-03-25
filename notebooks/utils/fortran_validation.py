"""
Fortran Code Validation Module

Provides functions to compare and validate generated Fortran files against reference files.
Used by both single_species_example.ipynb and multiple_species_example.ipynb.

**Normalization Strategy:**
The comparison applies selective normalization to ignore format variations that don't
affect code semantics, while preserving real code differences:

1. **Line Continuation Removal**
   - Removes Fortran line continuation markers ('&' or '+' at end of line)
   - Joins continued lines back together
   - Example: "a = b + &\n    c" → "a = b + c"

2. **Indentation Normalization (F77 Fixed Format)**
   - Strips leading whitespace from code lines (normalizes indentation variation)
   - Removes F77 labels (line numbers in columns 1-5)
   - Preserves comment lines as-is (starting with c, C, *, !)
   - Ensures consistency while respecting F77 column conventions

3. **Comment Whitespace Normalization**
   - Normalizes multiple spaces within comments to single spaces
   - Preserves comment content and structure
   - Examples: "c  Test" → "c Test", "c     Param  Data" → "c Param Data"

4. **Whitespace & Syntax Normalizations**
   - Whitespace around '=' (assignment statements)
   - Whitespace around mathematical operators: '+', '-', '*', '/'
   - Whitespace after commas in function calls
   - Fortran syntax: 'endif' → 'end if'
   - Fortran end statements: 'end subroutine name' → 'end'
   - Fortran operators: '.AND.', '.OR.', etc. → lowercase
   - Integer literals in assignments: '0' → '0.D0' (Fortran double precision)

Real code differences (different values, keywords, logic) are NOT ignored and will be reported.

Functions:
  - compare_fortran_files(): Direct comparison with selective normalization
  - validate_generated_fortran(): Batch validation of multiple files
  - print_validation_report(): Formatted validation report
"""

import re
from pathlib import Path
from difflib import unified_diff
from typing import Tuple, Dict, List
from decimal import Decimal, InvalidOperation


def compare_fortran_files(
    generated_file: Path,
    reference_file: Path,
    verbose: bool = False
) -> Tuple[bool, str]:
    """
    Compares two Fortran files DIRECTLY with selective normalization.

    This approach expects the macrofor library to generate code that exactly
    matches the reference. Minimal normalizations are applied only for
    format variations that don't affect code semantics:

    **Indentation Normalization (F77 Fixed Format):**
    - Strips leading whitespace from code lines (normalizes indentation variation)
    - Preserves comment lines as-is (starting with c, C, *, !)
    - Ensures consistency while respecting F77 column conventions

    **Whitespace & Syntax Normalizations:**
    - Whitespace around '=' (assignment statements only)
    - Whitespace around '/' (division operator)
    - Whitespace after commas in function calls
    - Fortran syntax: 'endif' → 'end if'
    - Fortran operators: '.AND.' → '.and.' (lowercase)
    - Fortran keywords: '.EQ.', '.NE.', '.LT.', '.LE.', '.GT.', '.GE.'
    - Integer literals in assignments: '0' → '0.D0' (Fortran double precision)

    Args:
        generated_file: Path to generated Fortran file
        reference_file: Path to reference Fortran file
        verbose: Show detailed diff output if files differ

    Returns:
        Tuple (is_identical: bool, error_message: str)
        - If identical: (True, "")
        - If different: (False, "error description")
    """
    if not generated_file.exists():
        return False, f"Generated file not found: {generated_file.name}"

    if not reference_file.exists():
        return False, f"Reference file not found: {reference_file.name}"

    # Read files with consistent encoding
    gen_content = generated_file.read_text(encoding='utf-8', errors='ignore')
    ref_content = reference_file.read_text(encoding='utf-8', errors='ignore')

    def remove_line_continuations(text: str) -> str:
        """
        Remove Fortran line continuations for both free- and fixed-format source.

        Rules implemented:
        - Free format: previous line ends with '&' or '+'. The marker is removed and the
          next line (optionally starting with '&' or '+') is concatenated without adding
          extra spaces.
        - Fixed format (F77): a line that has a non-blank continuation marker in column 6
          (index 5) continues the previous non-comment line. Columns 1-6 are removed and
          columns 7+ are concatenated to the previous line. Comment lines (c/C/*/!) are
          never treated as continuations.
        - Lines not marked as continuations start new records.
        """
        lines = text.split('\n')
        result: list[str] = []
        current: str | None = None

        for line in lines:
            stripped = line.lstrip()

            # Comment lines are never continuations
            is_comment = stripped.startswith(('c', 'C', '!', '*'))

            # Detect free-format continuation from previous line
            prev_continues_free = bool(
                current and current.rstrip().endswith(('&', '+'))
            )

            # Detect fixed-format continuation marker in THIS line (column 6)
            col6_marker = line[5] if len(line) > 5 else ' '
            is_fixed_cont_line = (not is_comment) and col6_marker not in (' ', '\t', '0')

            if prev_continues_free:
                # Remove trailing marker from previous line and optional marker from this line
                current = current.rstrip()[:-1]  # drop & or +
                # Remove leading continuation symbol if present on the continued line
                cont_payload = stripped
                if cont_payload.startswith(('&', '+')):
                    cont_payload = cont_payload[1:].lstrip()
                current = current + cont_payload
                continue

            if is_fixed_cont_line and current is not None:
                # Fixed format: drop cols 1-6, append cols 7+
                payload = line[6:] if len(line) > 6 else ''
                current = current + payload.lstrip()
                continue

            # If we reach here, this line starts a new record
            if current is not None:
                result.append(current)
            current = line

        if current is not None:
            result.append(current)

        return '\n'.join(result)

    def normalize_indentation(text: str) -> str:
        """
        Normalize Fortran F77 fixed-format indentation and remove labels.
        
        F77 format rules:
        - Columns 1-5: Line numbers/labels (usually empty)
        - Column 6: Continuation character (C, c, *, or !)
        - Columns 7-72: Code content
        - Columns 73-80: Comments/references (ignored)
        
        This function:
        1. Preserves comment lines (starting with c, C, *, !)
        2. Removes labels (columns 1-5) from code lines
        3. Strips leading whitespace from code lines (normalize indentation)
        4. Preserves code structure and content
        """
        lines = text.split('\n')
        normalized_lines = []
        
        for line in lines:
            # Check if it's a comment line (column 1 or after leading space)
            stripped = line.lstrip()
            if stripped.startswith(('!', 'c', 'C', '*')):
                # This is a comment - keep as-is
                normalized_lines.append(line)
            else:
                # Regular code line
                # Remove F77 labels (columns 1-5): digits followed by spaces
                # Pattern: line might start with up to 5 digits/spaces, then code
                line_no_label = re.sub(r'^\s{0,5}\d{1,5}\s*', '', line)
                
                # Strip leading whitespace to normalize indentation
                if line_no_label.strip():  # Non-empty line
                    normalized_lines.append(line_no_label.lstrip())
                else:  # Empty line
                    normalized_lines.append('')
        
        return '\n'.join(normalized_lines)

    def normalize_spacing(text: str) -> str:
        """Normalize whitespace and Fortran syntax variations."""
        def normalize_numeric_literals(line: str) -> str:
            """
            Normalize numerically equivalent literals to a canonical representation.

            Examples:
            - 1.D0 == 0.1D1
            - 0.1D-6 == 1D-7
            """

            number_pattern = re.compile(
                r'(?<![A-Za-z_])'
                r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[dDeE][-+]?\d+)?'
                r'(?![A-Za-z_])'
            )

            def repl(match: re.Match) -> str:
                token = match.group(0)
                # Leave pure integers untouched (array indices, labels, units)
                if re.fullmatch(r'[-+]?\d+', token):
                    return token

                normalized = token.replace('D', 'E').replace('d', 'e')
                try:
                    value = Decimal(normalized)
                except InvalidOperation:
                    return token

                # Canonical scientific notation using D-exponent
                if value == 0:
                    return '0.D0'

                s = f"{float(value):.15e}"  # e.g. 1.000000000000000e+00
                mantissa, exponent = s.split('e')
                mantissa = mantissa.rstrip('0').rstrip('.')
                if '.' not in mantissa:
                    mantissa += '.0'
                exp_int = int(exponent)
                return f"{mantissa}D{exp_int:+d}"

            return number_pattern.sub(repl, line)

        # First normalize indentation (F77 fixed format)
        text = normalize_indentation(text)
        
        lines = text.split('\n')
        normalized_lines = []

        for line in lines:
            # Strip ALL leading and trailing whitespace from already-normalized indentation
            # This ensures consistent spacing regardless of indentation differences
            line = line.strip()
            stripped = line
            # Check if it's a comment line: c, C, *, ! followed by space or nothing
            # (but not 'call', 'continue', etc.)
            is_comment = (
                stripped.startswith('!') or
                stripped.startswith('*') or
                (len(stripped) > 0 and stripped[0] in 'cC' and 
                 (len(stripped) == 1 or stripped[1] in ' \t'))
            )
            
            if is_comment:
                # Normalize whitespace in comment lines
                stripped = line.strip()
                # Empty comment line (just 'c' or 'c ' etc.) -> normalize to 'c'
                if len(stripped) <= 1 or (len(stripped) == 2 and stripped[1] == ' '):
                    line = stripped[0]  # Just 'c', '*', or '!'
                else:
                    # Remove extra spaces after comment character but keep content readable
                    line = re.sub(r'^(\s*[c*!])\s+', r'\1 ', line, flags=re.IGNORECASE)
                    # Compress multiple internal spaces to single spaces in comments
                    line = re.sub(r'(\s+)([^\s])', r' \2', line)
            else:
                # Code line - normalize whitespace
                
                # Normalize whitespace after commas (remove spaces)
                line = re.sub(r',\s+', ',', line)
                
                # Normalize whitespace around operators (=, +, -, *, /)
                line = re.sub(r'\s*=\s*', '=', line)
                line = re.sub(r'\s*\+\s*', '+', line)
                line = re.sub(r'\s*-\s*', '-', line)
                line = re.sub(r'\s*\*\s*', '*', line)
                line = re.sub(r'\s*/\s*', '/', line)

                # Normalize numeric zero literals in assignments: '=0' or '=0.0' → '=0.D0'
                # (spacing is already normalized, so match the compact form)
                line = re.sub(r'=0\.0(?!\d)', '=0.D0', line)
                line = re.sub(r'=0(?![\.\d])', '=0.D0', line)

                # Fortran 90 standard syntax
                line = re.sub(r'\bendif\b', 'end if', line, flags=re.IGNORECASE)
                
                # Normalize 'end subroutine <name>' to just 'end'
                line = re.sub(r'\bend\s+subroutine\s+\w+\s*$', 'end', line, flags=re.IGNORECASE)
                line = re.sub(r'\bend\s+function\s+\w+\s*$', 'end', line, flags=re.IGNORECASE)
                line = re.sub(r'\bend\s+program\s+\w+\s*$', 'end', line, flags=re.IGNORECASE)

                # Fortran logical operators (normalize to lowercase)
                line = re.sub(r'\.AND\.', '.and.', line, flags=re.IGNORECASE)
                line = re.sub(r'\.OR\.', '.or.', line, flags=re.IGNORECASE)
                line = re.sub(r'\.NOT\.', '.not.', line, flags=re.IGNORECASE)
                line = re.sub(r'\.EQ\.', '.eq.', line, flags=re.IGNORECASE)
                line = re.sub(r'\.NE\.', '.ne.', line, flags=re.IGNORECASE)
                line = re.sub(r'\.LT\.', '.lt.', line, flags=re.IGNORECASE)
                line = re.sub(r'\.LE\.', '.le.', line, flags=re.IGNORECASE)
                line = re.sub(r'\.GT\.', '.gt.', line, flags=re.IGNORECASE)
                line = re.sub(r'\.GE\.', '.ge.', line, flags=re.IGNORECASE)
                line = re.sub(r'\.TRUE\.', '.true.', line, flags=re.IGNORECASE)
                line = re.sub(r'\.FALSE\.', '.false.', line, flags=re.IGNORECASE)

                # Normalize numerically equivalent literals
                line = normalize_numeric_literals(line)

            normalized_lines.append(line)

        return '\n'.join(normalized_lines)

    # Remove line continuations FIRST (before any other normalization)
    gen_content = remove_line_continuations(gen_content)
    ref_content = remove_line_continuations(ref_content)
    
    gen_normalized = normalize_spacing(gen_content)
    ref_normalized = normalize_spacing(ref_content)

    # Direct comparison
    if gen_normalized == ref_normalized:
        return True, ""

    # Generate detailed error info if different
    gen_lines = gen_normalized.splitlines(keepends=True)
    ref_lines = ref_normalized.splitlines(keepends=True)

    if verbose:
        # Show unified diff
        diff = list(unified_diff(
            ref_lines,
            gen_lines,
            fromfile=f"Reference: {reference_file.name}",
            tofile=f"Generated: {generated_file.name}",
            lineterm=''
        ))
        diff_str = ''.join(diff[:100])  # First 100 lines
        return False, f"Files differ:\n{diff_str}"

    # Compact summary
    diff_count = sum(1 for gl, rl in zip(gen_lines, ref_lines) if gl != rl)
    line_diff = abs(len(gen_lines) - len(ref_lines))

    error_msg = f"Differences: {diff_count} lines differ"
    if line_diff > 0:
        error_msg += f", {line_diff} lines added/removed"

    return False, error_msg


def validate_generated_fortran(
    generated_dir: Path,
    reference_dir: Path,
    file_list: List[str] = None,
    verbose: bool = False
) -> Dict:
    """
    Validates all generated Fortran files against reference files.

    Args:
        generated_dir: Directory with generated Fortran files
        reference_dir: Directory with reference Fortran files
        file_list: List of specific files to check. If None, all *.f and *.inc files
        verbose: Detailed output (shows diffs)

    Returns:
        Dictionary with results:
        {
            'passed': [list of successful files],
            'failed': {filename: error_message},
            'missing': [generated files not in reference],
            'extra': [files in reference but not generated]
        }
    """
    results = {
        'passed': [],
        'failed': {},
        'missing': [],
        'extra': []
    }

    # Determine file list
    if file_list is None:
        gen_files = list(generated_dir.glob('*.f')) + list(generated_dir.glob('*.inc'))
        file_list = [f.name for f in gen_files]

    # Check each file
    for filename in file_list:
        gen_file = generated_dir / filename
        ref_file = reference_dir / filename

        if not gen_file.exists():
            results['missing'].append(filename)
            continue

        if not ref_file.exists():
            results['extra'].append(filename)
            continue

        is_identical, error_msg = compare_fortran_files(gen_file, ref_file, verbose)

        if is_identical:
            results['passed'].append(filename)
        else:
            results['failed'][filename] = error_msg

    return results


def print_validation_report(results: Dict, test_name: str = "Validation") -> bool:
    """
    Prints a formatted validation report.

    Args:
        results: Dictionary from validate_generated_fortran()
        test_name: Name of the test (for display)

    Returns:
        True if all files passed, False otherwise
    """
    total = len(results['passed']) + len(results['failed']) + len(results['missing'])

    print("\n" + "=" * 70)
    print(f"{test_name.upper()} - REPORT")
    print("=" * 70)

    # Passed tests
    if results['passed']:
        print(f"\n✓ PASSED ({len(results['passed'])}/{total}):")
        for fname in sorted(results['passed']):
            print(f"  ✓ {fname}")

    # Failed tests
    if results['failed']:
        print(f"\n✗ FAILED ({len(results['failed'])}/{total}):")
        for fname, msg in sorted(results['failed'].items()):
            print(f"  ✗ {fname}")
            if msg:
                first_line = msg.split('\n')[0]
                print(f"    → {first_line}")

    # Missing files
    if results['missing']:
        print(f"\n⚠ MISSING ({len(results['missing'])}):")
        for fname in sorted(results['missing']):
            print(f"  ⚠ {fname}")

    # Extra files
    if results['extra']:
        print(f"\n⊕ EXTRA ({len(results['extra'])}):")
        for fname in sorted(results['extra']):
            print(f"  ⊕ {fname}")

    # Summary
    print("\n" + "-" * 70)
    success_rate = (len(results['passed']) / total * 100) if total > 0 else 0
    print(f"TOTAL: {len(results['passed'])}/{total} passed ({success_rate:.1f}%)")
    print("=" * 70)

    return len(results['failed']) == 0 and len(results['missing']) == 0
