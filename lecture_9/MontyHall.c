/*
Monty Hall Simulator: A program to simulate the Monty Hall problem using A, B, C doors.
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h> // Included to handle case conversion for user input

int main()
{
    // Where is the grand prize?
    srand(time(NULL));

    int prize = rand() % 3;
    int notprize1 = (prize + 1) % 3;
    int notprize2 = (prize + 2) % 3;

    // Ask contestant to pick a door
    printf("\n");
    printf(" =-=-=-=-=-=-=-=-=-=-=-=- \n");
    printf(" ** Monty Hall Simulator ** \n");
    printf(" =-=-=-=-=-=-=-=-=-=-=-=- \n");
    printf("\n");

    char userChar;
    int pick = -1;

    printf(" Pick a door (A, B, or C): ");
    scanf(" %c", &userChar);

    userChar = toupper(userChar);

    // map characters to integers (A=0, B=1, C=2)
    if (userChar == 'A') pick = 0;
    else if (userChar == 'B') pick = 1;
    else if (userChar == 'C') pick = 2;

    // validate input
    while(pick < 0 || pick > 2)
    {
        printf(" Invalid selection. Please pick A, B, or C: ");
        scanf(" %c", &userChar);
        userChar = toupper(userChar);
        
        if (userChar == 'A') pick = 0;
        else if (userChar == 'B') pick = 1;
        else if (userChar == 'C') pick = 2;
        else pick = -1;
    }

    printf("\n You entered: %c\n", userChar);

    // Tell contestant about another door
    int other;       // The door Monty opens
    int other_other; // The remaining closed door 
    
    printf("\n Interesting choice ...\n");

    // Logic to decide which door to open
    if (pick == prize)
    {
        int ss = rand() % 2;
        if (ss == 0)
        {
            other = notprize1;
            other_other = notprize2;
        }
        else
        {
            other = notprize2;
            other_other = notprize1;
        }
    }
    else
    {
        // If contestant didn't pick the prize, Monty opens the only other loser
        other_other = prize; 
        
        if (pick == notprize1)
        {
            other = notprize2;
        }
        else
        {
            other = notprize1;
        }
    }

    // Convert 'other' int back to char
    char otherChar = 'A' + other;
    char switchChar = 'A' + other_other;

    printf("\n I can tell you for sure that the prize is not behind Door: %c\n", otherChar);

    // Ask if they want to change
    int change = -1;
    while(change != 0 && change != 1)
    {
        printf("\n Stay with Door %c (press 0) or switch to Door %c (press 1): ", userChar, switchChar);
        scanf("%d", &change);
    }

    // Determine final pick
    int final_pick;
    char finalChar;

    if (change == 0)
    {
        final_pick = pick;
        finalChar = userChar;
        printf("\n You stayed with Door %c\n", finalChar);
    }
    else
    {
        final_pick = other_other;
        finalChar = switchChar;
        printf("\n You switched to Door %c\n", finalChar);
    }

    // Check answer
    if (final_pick == prize)
    {
        printf("\n *** WINNER ***\n\n");
    }
    else
    {
        printf("\n --- LOSER ---\n\n");
    }

    printf(" The prize was behind Door %c\n\n", 'A' + prize);

    return 0;
}