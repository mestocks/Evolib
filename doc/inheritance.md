          msFormat       FastaFormat
              ^          ^         ^
              |          |         |
           IoPopulationData     DnaPopulationData ::- FastaSequence <-- DnaSequence
                  ^                     ^         ::- FastaSequence
                  |                     |         ::- ...
              IoStats               DnaStats
