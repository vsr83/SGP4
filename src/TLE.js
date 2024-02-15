const TleFieldType = {
    FIELD_NUMBER : 1,
    FIELD_NUMBER_IMPLIED : 2,
    FIELD_NUMBER_MAPPED : 3,
    FIELD_STRING : 4
};

const tleFormat = [
    {name : "title",         fieldType : TleFieldType.FIELD_STRING, row : 0, startCol :  0, endCol : 23},
    {name : "catalogNumber", fieldType : TleFieldType.FIELD_NUMBER, row : 1, startCol :  2, endCol :  6},
    {name : "classification",fieldType : TleFieldType.FIELD_STRING, row : 1, startCol :  7, endCol :  7},
    {name : "intLaunchYear", fieldType : TleFieldType.FIELD_STRING, row : 1, startCol :  9, endCol : 10},
    {name : "intLaunchNum",  fieldType : TleFieldType.FIELD_STRING, row : 1, startCol : 11, endCol : 13},
    {name : "intLaunchPiece",fieldType : TleFieldType.FIELD_STRING, row : 1, startCol : 14, endCol : 16},
    {name : "epochYear",     fieldType : TleFieldType.FIELD_NUMBER_MAPPED, row : 1, startCol :18, endCol :19,
     callback : function(suffix) 
     {
        if (suffix > 56) { 
            return 1900 + suffix;
        }
        else 
        {
            return 2000 + suffix;
        }
    }},
    {name : "epochFracDay",     fieldType : TleFieldType.FIELD_NUMBER, row : 1, startCol : 20, endCol : 31},
    {name : "meanMotionDer",    fieldType : TleFieldType.FIELD_NUMBER, row : 1, startCol : 33, endCol : 42},
    {name : "meanMotionDer2",   fieldType : TleFieldType.FIELD_NUMBER_IMPLIED, row : 1, startCol : 44, endCol : 51},
    {name : "dragTerm",         fieldType : TleFieldType.FIELD_NUMBER_IMPLIED, row : 1, startCol : 53, endCol : 60},
    {name : "ephemerisType",    fieldType : TleFieldType.FIELD_NUMBER, row : 1, startCol : 62, endCol : 62},
    {name : "elementSetNo",     fieldType : TleFieldType.FIELD_NUMBER, row : 1, startCol : 64, endCol : 67},
    {name : "checkSum1",        fieldType : TleFieldType.FIELD_NUMBER, row : 1, startCol : 68, endCol : 68},
    {name : "catalogNumber2",   fieldType : TleFieldType.FIELD_NUMBER, row : 2, startCol :  2, endCol :  6},
    {name : "inclination",      fieldType : TleFieldType.FIELD_NUMBER, row : 2, startCol :  8, endCol : 15},
    {name : "raAscNode",        fieldType : TleFieldType.FIELD_NUMBER, row : 2, startCol : 17, endCol : 24},
    {name : "eccentricity",     fieldType : TleFieldType.FIELD_NUMBER_IMPLIED, row : 2, startCol : 26, endCol : 32},
    {name : "argPerigee",       fieldType : TleFieldType.FIELD_NUMBER, row : 2, startCol : 34, endCol : 41},
    {name : "meanAnomaly",      fieldType : TleFieldType.FIELD_NUMBER, row : 2, startCol : 43, endCol : 50},
    {name : "meanMotion",       fieldType : TleFieldType.FIELD_NUMBER, row : 2, startCol : 52, endCol : 62},
    {name : "revNoAtEpoch",     fieldType : TleFieldType.FIELD_NUMBER, row : 2, startCol : 63, endCol : 67},
    {name : "checkSum2",        fieldType : TleFieldType.FIELD_NUMBER, row : 2, startCol : 68, endCol : 68}
];

/**
 * Create from lines.
 * 
 * @param {string[]} lines
 *      The three lines of the TLE as strings.
 * @returns {Tle} The object constructed from the lines.
 */
export function tleFromLines(lines) 
{
    const tle = {};

    /**
     * Fill from a JSON.
     * 
     * @param {ITle} json 
     *      The JSON.
     */
    function fillFromJson(json)
    {
        tle.title = json.title;
        tle.catalogNumber = json.catalogNumber;
        tle.classification = json.classification;
        tle.intLaunchYear = json.intLaunchYear;
        tle.intLaunchNum = json.intLaunchNum;
        tle.intLaunchPiece = json.intLaunchPiece;
        tle.epochYear = json.epochYear;
        tle.epochFracDay = json.epochFracDay;
        tle.meanMotionDer = json.meanMotionDer;
        tle.meanMotionDer2 = json.meanMotionDer2;
        tle.dragTerm = json.dragTerm;
        tle.ephemerisType = json.ephemerisType;
        tle.elementSetNo = json.elementSetNo;
        tle.checkSum1 = json.checkSum1;
        tle.catalogNumber2 = json.catalogNumber2;
        tle.inclination = json.inclination;
        tle.raAscNode = json.raAscNode;
        tle.eccentricity = json.eccentricity;
        tle.argPerigee = json.argPerigee;
        tle.meanAnomaly = json.meanAnomaly;
        tle.meanMotion = json.meanMotion;
        tle.revNoAtEpoch = json.revNoAtEpoch;
        tle.checkSum2 = json.checkSum2;
    }

    const json = {};

    for (let indField = 0; indField < tleFormat.length; indField++)
    {
        const field = tleFormat[indField];
        const str = lines[field.row].substring(field.startCol, field.endCol + 1);

        switch(field.fieldType)
        {
        case TleFieldType.FIELD_STRING:
            json[field.name] = str;
            break;
        case TleFieldType.FIELD_NUMBER:
            json[field.name] = Number(str);
            break;
        case TleFieldType.FIELD_NUMBER_IMPLIED:
            if (str[0] == '-')
            {
                json[field.name] = -Number("0." + str.substring(1).replace("-", "e-").replace("+", "e").replace(" ", ""));
            }
            else 
            {
                json[field.name] = Number("0." + str.replace("-", "e-").replace("+", "e").replace(" ", ""));
            }
            break;
        case TleFieldType.FIELD_NUMBER_MAPPED:
            if (field.callback)
            {
                json[field.name] = field.callback(Number(str));
            }
            break;
        }
    }
    fillFromJson(json);

    tle.checkSumValid = (tleLineChecksum(lines[1]) == tle.checkSum1) &&
                        (tleLineChecksum(lines[2]) == tle.checkSum2);

    tle.jtUt1Epoch = tleParseEpoch(tle.epochYear, tle.epochFracDay);

    return tle;
}

/**
     * Compute checksum for TLE line.
     * 
     * @param {string} line 
     *      TLE line
     * @returns {number} Modulo-10 checksum.
     */
function tleLineChecksum(lineIn)
{
    let checksum = 0; 
    const line = lineIn.substring(0, 68);

    for (let ind = 0; ind < line.length; ind++)
    {
        const char = line[ind];

        if ('0123456789'.indexOf(char) > -1)
        {
            checksum += parseInt(char);
        }
        else if (char == '-')
        {
            checksum++;
        }
    }
    return checksum % 10;
}

/**
 * Compute Julian time of the epoch.
 * 
 * @param {number} year 
 *      The year (all four digits).
 * @param {number} days 
 *      Fractional days of the year.
 * @returns {number} The Julian time of the epoch.
 */
function tleParseEpoch(year, days)
{
    const isLeap = (year % 4) == 0;
    const monthLengths = [31, isLeap ? 29 : 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    const dayOfYear = Math.floor(days);

    let month = 0;
    let sumDays = 0;
    for (month = 0; month < 12; month++)
    {
        const monthLength = monthLengths[month];

        if (dayOfYear - sumDays - monthLength < 0)
        {
            break;
        }
        sumDays += monthLength;
    }
    const dayOfMonth = dayOfYear - sumDays;
    const fracDay = days - dayOfYear;
    const hour = Math.floor(fracDay * 24.0);
    const minute = Math.floor((fracDay - hour / 24.0) * 1440.0);
    const second = (fracDay - hour / 24.0 - minute / 1440.0) * 86400.0;

    return timeJulianYmdhms(year, month + 1, dayOfMonth, hour, minute, second);
}

/**
 * Compute Julian date for given calendar date.
 * 
 * @param {number} year 
 *      Year as an integer.
 * @param {number} month 
 *      Month (1-12).
 * @param {number} mday 
 *      Day of the month (1-31).
 * @returns {number} Julian date.
 */
function dateJulianYmd(year, month, mday)
{
    if (month < 3)
    {
        year--;
        month += 12;
    }

    const A = Math.floor(year / 100.0);
    const B = Math.floor(A / 4.0);
    const C = Math.floor(2.0 - A + B);
    const E = Math.floor(365.25 * (year + 4716.0));
    const F = Math.floor(30.6001 * (month + 1.0));

    return C + mday + E + F - 1524.5;    
}

/**
 * Compute Julian time.
 * 
 * @param {number} year 
 *      Year as an integer.
 * @param {number} month 
 *      Month (1-12) integer.
 * @param {number} mday 
 *      Day of the month (1-31) integer.
 * @param {number} hour 
 *      Hour (0-23) integer.
 * @param {number} minute
 *      Minute (0-59) integer. 
 * @param {number} second 
 *      Second (0-60) floating point.
 * @returns {number} An object with JD and JT for Julian date and time.
 */
export function timeJulianYmdhms(year, month, mday, hour, minute, second)
{
    const JD = dateJulianYmd(year, month, mday);
    const JT = JD + hour / 24.0 + minute/(24.0 * 60.0) + second/(24.0 * 60.0 * 60.0);

    return JT;
}