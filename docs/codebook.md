# Data codebook

This file maps the sanitised R column names used in the scripts to the variable descriptions reported in the manuscript.

The manuscript states that the study dataset contains 46 columns and 866 daily rows, covering 2020-01-30 to 2022-06-13. It also describes the variables listed below.

| R column name | Paper label | Description |
|---|---|---|
| `Date` | Date | Date observations were recorded. |
| `Excess.mortality` | Excess mortality | Percentage difference between reported and projected deaths. |
| `Schools` | Schools | Whether schools were open, partially open, or closed. |
| `Face.masks` | Face masks | Mask mandates during the pandemic: optional, mandatory, or no mandate. |
| `Lockdown.severity` | Lockdown severity | Severity of lockdown: severe, moderate, weak, or limited measures/social distancing only. |
| `Majority.COVID.19.variant` | Majority COVID-19 variant | The majority COVID-19 variant. |
| `Flights.7.day.moving.average` | Flights 7-day moving average | 7-day moving average of flights. |
| `OpenTable.restaurant.bookings.London.index` | OpenTable restaurant | OpenTable index on restaurant bookings in London. |
| `Google.homeworking.Greater.London.mobility.index` | Google homeworking | Google index on homeworking in Greater London. |
| `Google.workplace.Greater.London.mobility.index` | Google workplace | Google index on workplace activity in Greater London. |
| `Apple.walking.London.mobility.index` | Apple walking | Apple index on walking activity in London. |
| `Google.parks.Greater.London.mobility.index` | Google parks | Google index on park visits in Greater London. |
| `Google.retail.recreation.Greater.London.mobility.index` | Google retail & recreation | Google index on retail and recreation in Greater London. |
| `Google.grocery.pharmacy.Greater.London.mobility.index` | Google grocery & pharmacy | Google index on grocery and pharmacy in Greater London. |
| `Google.transit.stations.mobility.index` | Google transit stations | Google index on transit station activity. |
| `TfL.Tube.mobility.index` | TfL Tube | TfL index on tube activity. |
| `TfL.Bus.mobility.index` | TfL Bus | TfL index on bus activity. |
| `Citymapper.journeys.mobility.index` | Citymapper journeys | Citymapper index on journey activity. |
| `Season` | Season | Winter, autumn, summer, and spring. |
| `PCR.tests` | PCR tests | Number of tests on date. |
| `PCR.tests.capacity` | PCR tests capacity | Capacity level on date. |
| `Antibody.tests` | Antibody tests | Number of tests on date. |
| `Antibody.tests.capacity` | Antibody tests capacity | Capacity level on date. |
| `Pillar.1.NHS.and.UKHSA.capacity` | Pillar 1 capacity | NHS/UKHSA capacity on date. |
| `Pillar.2.UK.Government.capacity` | Pillar 2 capacity | UK Government capacity on date. |
| `Pillar.3.Antibody.capacity` | Pillar 3 capacity | Antibody testing capacity on date. |
| `Pillar.4.Surveillance.capacity` | Pillar 4 capacity | Surveillance testing capacity on date. |
| `Pillar.1.NHS.and.UKHSA.tests` | Pillar 1 tests | NHS/UKHSA tests on date. |
| `Pillar.2.UK.Government.tests` | Pillar 2 tests | UK Government tests on date. |
| `Pillar.3.Antibody.tests` | Pillar 3 tests | Antibody tests on date. |
| `Pillar.4.Surveillance.tests` | Pillar 4 tests | Surveillance tests on date. |
| `Tests.across.all.4.Pillars` | Tests across all 4 Pillars | Total tests on date. |
| `New.cases` | New cases | People testing positive for COVID-19 on date. |
| `New.infections` | New infections | New infections on date. |
| `Reinfections` | Reinfections | New reinfections on date. |
| `Hospital.admissions` | Hospital admissions | Patients admitted to hospital with COVID-19. |
| `Patients.in.hospital` | Patients in hospital | Patients in hospital with COVID-19. |
| `Patients.in.MVBs` | COVID-19 patients in MVBs | Patients in Mechanical Ventilator Beds (MVBs) with COVID-19. |
| `Vaccinations.total` | Vaccinations (total) | Total vaccines administered on date. |
| `Vaccinations.1st.dose` | Vaccinations (1st dose) | First dose vaccines administered on date. |
| `Vaccinations.2nd.dose` | Vaccinations (2nd dose) | Second dose vaccines administered on date. |
| `Vaccinations.3rd.dose` | Vaccinations (3rd dose) | Third dose vaccines administered on date. |
| `First.dose.uptake` | 1st dose uptake | Reported first dose uptake. |
| `Second.dose.uptake` | 2nd dose uptake | Reported second dose uptake. |
| `Third.dose.uptake` | 3rd dose uptake | Reported third dose uptake. |
| `COVID.19.deaths.on.certificate` | COVID-19 deaths on certificate | Daily deaths with COVID-19 on certificate by date of death. |

## Columns removed in preprocessing

The uploaded read script drops three columns before analysis:

- `New.cases.specimen.date.7.day.change.`
- `New.cases.specimen.date.7.day...change.`
- `Deaths.within.28.days.of.COVID.19` 
